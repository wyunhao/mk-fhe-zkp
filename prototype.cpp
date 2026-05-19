// Multi-key RLWE Vector Encryption (BGV) — prototype of one layer of multi-key multiplication.
//
// Implements Definition 5.2 (RLWE-based Vector Encryption Construction):
//   Setup, E, D, mkMul, mkDec — over BGV with ternary HW-bounded secrets.
//
// Parameters (per spec):
//   N = 8192, p = 7*2^20 + 1, log q ~ 130 bits (three small primes), manual mod-switch to log q' ~ 44 bits,
//   secret key components are ternary with fixed Hamming weight h = 32.
//
// Notation in code:
//   * `s` for party 1 (scalar enc):  single secret polynomial in R_q.
//   * `s_vec[i]`, i in [ell'] for party 2 (vector enc): ell' = ell + tau parallel secrets sharing a common `a`.
//   * `T`: sparsification matrix in R_p^{tau x ell}, used to bind a redundancy check on the encoded vector.
//
// Built on top of microsoft/SEAL's low-level util layer (NTT/dyadic mul/RNS arithmetic).

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "seal/memorymanager.h"
#include "seal/modulus.h"
#include "seal/util/ntt.h"
#include "seal/util/numth.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/uintarithsmallmod.h"

using namespace seal;
using namespace seal::util;
using std::size_t;
using std::uint64_t;
using std::vector;

namespace mk {

// -------------------------- PRNG --------------------------
class Prng {
public:
    Prng() : rng_(std::random_device{}()) {}
    explicit Prng(uint64_t seed) : rng_(seed) {}
    uint64_t u64() { return rng_(); }
    uint64_t uniform_below(uint64_t bound) { return bound ? (rng_() % bound) : 0; }
    // Centered binomial distribution: sum(eta bits) - sum(eta bits). sigma = sqrt(eta/2).
    int cbd(int eta) {
        int s = 0;
        for (int i = 0; i < eta; i++) s += int(rng_() & 1ULL);
        for (int i = 0; i < eta; i++) s -= int(rng_() & 1ULL);
        return s;
    }
private:
    std::mt19937_64 rng_;
};

// -------------------------- Params --------------------------
struct Params {
    size_t N = 0;
    int log_N = 0;
    Modulus p;                                       // plaintext modulus
    vector<Modulus> q_chain;                         // q_0, q_1, ..., q_{L-1}; mod-switch drops the back
    vector<std::unique_ptr<NTTTables>> q_ntt;
    std::unique_ptr<NTTTables> p_ntt;                // for plaintext-side R_p ops (p is NTT-friendly)
    size_t ell = 0, tau = 0, ell_prime = 0;
    size_t hw = 0;                                   // Hamming weight for ternary secret
    vector<uint64_t> p_mod_q;                        // p reduced modulo each q_i
    int sigma_eta = 21;                              // CBD parameter for noise; sigma ~= sqrt(21/2) ~= 3.24
};

static void params_init(Params& params,
                        uint64_t plain_modulus_val,
                        const vector<int>& q_bit_sizes,
                        size_t N, size_t ell, size_t tau, size_t hw)
{
    params.N = N;
    params.log_N = 0;
    while ((size_t(1) << params.log_N) < N) params.log_N++;
    if ((size_t(1) << params.log_N) != N) throw std::runtime_error("N must be power of 2");

    params.p = Modulus(plain_modulus_val);
    params.ell = ell;
    params.tau = tau;
    params.ell_prime = ell + tau;
    params.hw = hw;

    // Generate q primes, each ~ bit_size and ≡ 1 (mod 2N) so NTT exists.
    // Disambiguate duplicates if same bit_size repeats.
    params.q_chain.clear();
    for (int bs : q_bit_sizes) {
        auto candidates = get_primes(2 * uint64_t(N), bs, 8);
        Modulus chosen(0);
        for (auto& c : candidates) {
            bool dup = false;
            for (auto& q : params.q_chain) if (q.value() == c.value()) { dup = true; break; }
            if (!dup) { chosen = c; break; }
        }
        if (chosen.value() == 0) throw std::runtime_error("could not allocate q prime of given bit size");
        params.q_chain.push_back(chosen);
    }

    params.q_ntt.clear();
    params.q_ntt.reserve(params.q_chain.size());
    for (auto& q : params.q_chain) {
        params.q_ntt.emplace_back(std::make_unique<NTTTables>(params.log_N, q));
    }

    // plaintext NTT (we need p ≡ 1 mod 2N; for our p = 7*2^20+1 with 2N=16384, this holds since 2^14 | 7*2^20).
    params.p_ntt = std::make_unique<NTTTables>(params.log_N, params.p);

    params.p_mod_q.resize(params.q_chain.size());
    for (size_t i = 0; i < params.q_chain.size(); i++) {
        params.p_mod_q[i] = barrett_reduce_64(plain_modulus_val, params.q_chain[i]);
    }
}

// -------------------------- RNS polynomial --------------------------
// Layout: prime-major. data has size num_primes * N. Coefficient i of prime k at data[k*N + i].
// All polys are stored in NTT form unless explicitly noted.
struct RnsPoly {
    vector<uint64_t> data;
    size_t num_primes = 0;
    size_t N = 0;
};

static inline uint64_t* row(RnsPoly& p, size_t k) { return p.data.data() + k * p.N; }
static inline const uint64_t* row(const RnsPoly& p, size_t k) { return p.data.data() + k * p.N; }

static RnsPoly make_zero(const Params& params, size_t num_primes) {
    RnsPoly r;
    r.num_primes = num_primes;
    r.N = params.N;
    r.data.assign(num_primes * params.N, 0ULL);
    return r;
}

// Sample uniform R_q directly in NTT form (NTT preserves the uniform distribution).
static RnsPoly sample_uniform(const Params& params, size_t num_primes, Prng& prng) {
    RnsPoly r = make_zero(params, num_primes);
    for (size_t k = 0; k < num_primes; k++) {
        uint64_t q = params.q_chain[k].value();
        uint64_t* dst = row(r, k);
        for (size_t i = 0; i < params.N; i++) dst[i] = prng.u64() % q;
    }
    return r;
}

// Sample HW=h ternary in {-1,0,1}^N — coefficient form, return signed ints.
static vector<int> sample_ternary_hw(size_t N, size_t hw, Prng& prng) {
    if (hw > N) throw std::runtime_error("hw > N");
    vector<int> coeffs(N, 0);
    vector<size_t> positions(N);
    for (size_t i = 0; i < N; i++) positions[i] = i;
    for (size_t i = 0; i < hw; i++) {
        size_t j = i + prng.uniform_below(N - i);
        std::swap(positions[i], positions[j]);
    }
    for (size_t i = 0; i < hw; i++) {
        coeffs[positions[i]] = (prng.u64() & 1ULL) ? 1 : -1;
    }
    return coeffs;
}

// Convert (small) signed coeffs to RNS-NTT form at given level.
static RnsPoly small_signed_to_rns_ntt(const vector<int>& coeffs, const Params& params, size_t num_primes) {
    if (coeffs.size() != params.N) throw std::runtime_error("size mismatch");
    RnsPoly r = make_zero(params, num_primes);
    for (size_t k = 0; k < num_primes; k++) {
        uint64_t q = params.q_chain[k].value();
        uint64_t* dst = row(r, k);
        for (size_t i = 0; i < params.N; i++) {
            int v = coeffs[i];
            dst[i] = (v >= 0) ? uint64_t(v) : (q - uint64_t(-v));
        }
        ntt_negacyclic_harvey(dst, *params.q_ntt[k]);
    }
    return r;
}

// Lift plaintext (coeffs in [0,p)) into R_q in NTT form for each prime.
static RnsPoly plain_to_rns_ntt(const vector<uint64_t>& coeffs_p, const Params& params, size_t num_primes) {
    if (coeffs_p.size() != params.N) throw std::runtime_error("size mismatch");
    RnsPoly r = make_zero(params, num_primes);
    for (size_t k = 0; k < num_primes; k++) {
        uint64_t* dst = row(r, k);
        // coeffs_p[i] < p < q_k, copy directly
        std::memcpy(dst, coeffs_p.data(), params.N * sizeof(uint64_t));
        ntt_negacyclic_harvey(dst, *params.q_ntt[k]);
    }
    return r;
}

// Sample CBD noise coefficients.
static vector<int> sample_noise_cbd(size_t N, int eta, Prng& prng) {
    vector<int> coeffs(N);
    for (size_t i = 0; i < N; i++) coeffs[i] = prng.cbd(eta);
    return coeffs;
}

// -------------------------- Poly arithmetic (NTT form, dyadic) --------------------------
static RnsPoly add(const RnsPoly& a, const RnsPoly& b, const Params& params) {
    if (a.num_primes != b.num_primes || a.N != b.N) throw std::runtime_error("level mismatch in add");
    RnsPoly r = make_zero(params, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++) {
        add_poly_coeffmod(row(a, k), row(b, k), a.N, params.q_chain[k], row(r, k));
    }
    return r;
}

static RnsPoly sub(const RnsPoly& a, const RnsPoly& b, const Params& params) {
    if (a.num_primes != b.num_primes || a.N != b.N) throw std::runtime_error("level mismatch in sub");
    RnsPoly r = make_zero(params, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++) {
        sub_poly_coeffmod(row(a, k), row(b, k), a.N, params.q_chain[k], row(r, k));
    }
    return r;
}

static RnsPoly mul(const RnsPoly& a, const RnsPoly& b, const Params& params) {
    if (a.num_primes != b.num_primes || a.N != b.N) throw std::runtime_error("level mismatch in mul");
    RnsPoly r = make_zero(params, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++) {
        dyadic_product_coeffmod(row(a, k), row(b, k), a.N, params.q_chain[k], row(r, k));
    }
    return r;
}

static RnsPoly mul_by_p(const RnsPoly& a, const Params& params) {
    RnsPoly r = make_zero(params, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++) {
        multiply_poly_scalar_coeffmod(row(a, k), a.N, params.p_mod_q[k], params.q_chain[k], row(r, k));
    }
    return r;
}

// Truncate a full-level secret to a lower level (drop tail primes).
static RnsPoly truncate_to_level(const RnsPoly& s, const Params& params, size_t target_primes) {
    if (target_primes > s.num_primes) throw std::runtime_error("truncate up not allowed");
    RnsPoly r;
    r.N = params.N;
    r.num_primes = target_primes;
    r.data.assign(s.data.begin(), s.data.begin() + target_primes * params.N);
    return r;
}

// -------------------------- BGV mod-switch (drop last prime) --------------------------
// Operates on one R_q polynomial c (in NTT form, with .num_primes = L).
// After: c is at level L-1. Mirrors SEAL's mod_t_and_divide_q_last_ntt_inplace.
// Effect on plaintext mod p: c' ≡ c · q_last^{-1} (mod p)
static void mod_switch_drop_last_poly(RnsPoly& c, const Params& params) {
    if (c.num_primes <= 1) throw std::runtime_error("cannot drop last prime: only one left");
    const size_t L = c.num_primes;
    const Modulus& q_last = params.q_chain[L - 1];
    const uint64_t q_last_val = q_last.value();

    // Step 1: c_last in coeff form
    vector<uint64_t> c_last(params.N);
    std::memcpy(c_last.data(), row(c, L - 1), params.N * sizeof(uint64_t));
    inverse_ntt_negacyclic_harvey(c_last.data(), *params.q_ntt[L - 1]);

    // Step 2: neg_c_last_mod_p = (-c_last mod p), then multiplied by q_last^{-1} mod p (skip if =1)
    vector<uint64_t> ncl_mod_p(params.N);
    modulo_poly_coeffs(c_last.data(), params.N, params.p, ncl_mod_p.data());
    negate_poly_coeffmod(ncl_mod_p.data(), params.N, params.p, ncl_mod_p.data());

    uint64_t inv_q_last_mod_p = 1;
    if (!try_invert_uint_mod(q_last_val, params.p, inv_q_last_mod_p))
        throw std::runtime_error("q_last not invertible mod p");
    if (inv_q_last_mod_p != 1) {
        multiply_poly_scalar_coeffmod(ncl_mod_p.data(), params.N, inv_q_last_mod_p, params.p, ncl_mod_p.data());
    }

    vector<uint64_t> delta(params.N);
    for (size_t k = 0; k < L - 1; k++) {
        const Modulus& q_k = params.q_chain[k];
        const uint64_t q_k_val = q_k.value();

        // delta = ncl_mod_p mod q_k (coeff form)
        modulo_poly_coeffs(ncl_mod_p.data(), params.N, q_k, delta.data());

        // delta *= q_last (mod q_k)
        uint64_t qlast_mod_qk = barrett_reduce_64(q_last_val, q_k);
        multiply_poly_scalar_coeffmod(delta.data(), params.N, qlast_mod_qk, q_k, delta.data());

        // delta += c_last mod q_k  (coeff form)
        for (size_t i = 0; i < params.N; i++) {
            uint64_t cl_mod_qk = barrett_reduce_64(c_last[i], q_k);
            delta[i] = add_uint_mod(delta[i], cl_mod_qk, q_k);
        }

        // NTT delta; subtract from c_k (in NTT form)
        ntt_negacyclic_harvey(delta.data(), *params.q_ntt[k]);
        sub_poly_coeffmod(row(c, k), delta.data(), params.N, q_k, row(c, k));

        // c_k *= q_last^{-1} (mod q_k)
        uint64_t inv_qlast_mod_qk = 0;
        if (!try_invert_uint_mod(q_last_val, q_k, inv_qlast_mod_qk))
            throw std::runtime_error("q_last not invertible mod q_k");
        multiply_poly_scalar_coeffmod(row(c, k), params.N, inv_qlast_mod_qk, q_k, row(c, k));
    }

    // Drop last prime
    c.data.resize((L - 1) * params.N);
    c.num_primes = L - 1;
}

// -------------------------- Plaintext arithmetic in R_p --------------------------
static vector<uint64_t> plain_mul(const vector<uint64_t>& a, const vector<uint64_t>& b, const Params& params) {
    vector<uint64_t> an(a), bn(b), c(params.N);
    ntt_negacyclic_harvey(an.data(), *params.p_ntt);
    ntt_negacyclic_harvey(bn.data(), *params.p_ntt);
    dyadic_product_coeffmod(an.data(), bn.data(), params.N, params.p, c.data());
    inverse_ntt_negacyclic_harvey(c.data(), *params.p_ntt);
    return c;
}

static vector<uint64_t> plain_add(const vector<uint64_t>& a, const vector<uint64_t>& b, const Params& params) {
    vector<uint64_t> c(params.N);
    add_poly_coeffmod(a.data(), b.data(), params.N, params.p, c.data());
    return c;
}

// -------------------------- Secret keys --------------------------
struct SecretKeyScalar {
    vector<int> s_coeffs;     // signed ternary, |s_coeffs[i]| in {0,1}
    RnsPoly s_ntt;            // full-level NTT-RNS
};

struct SecretKeyVector {
    // ell' parallel secrets, each ternary HW=hw, all sharing a common a (the a lives inside the ciphertext)
    vector<vector<int>> s_coeffs;
    vector<RnsPoly> s_ntt;
    // Sparsification matrix T in R_p^{tau x ell}. T_coeffs[r][c] is a poly's coeffs in [0,p).
    vector<vector<vector<uint64_t>>> T_coeffs;
};

static SecretKeyScalar keygen_scalar(const Params& params, Prng& prng) {
    SecretKeyScalar sk;
    sk.s_coeffs = sample_ternary_hw(params.N, params.hw, prng);
    sk.s_ntt = small_signed_to_rns_ntt(sk.s_coeffs, params, params.q_chain.size());
    return sk;
}

static SecretKeyVector keygen_vector(const Params& params, Prng& prng) {
    SecretKeyVector sk;
    sk.s_coeffs.resize(params.ell_prime);
    sk.s_ntt.reserve(params.ell_prime);
    for (size_t i = 0; i < params.ell_prime; i++) {
        sk.s_coeffs[i] = sample_ternary_hw(params.N, params.hw, prng);
        sk.s_ntt.emplace_back(small_signed_to_rns_ntt(sk.s_coeffs[i], params, params.q_chain.size()));
    }
    sk.T_coeffs.assign(params.tau, vector<vector<uint64_t>>(params.ell, vector<uint64_t>(params.N, 0)));
    for (size_t r = 0; r < params.tau; r++) {
        for (size_t c = 0; c < params.ell; c++) {
            for (size_t i = 0; i < params.N; i++) {
                sk.T_coeffs[r][c][i] = prng.u64() % params.p.value();
            }
        }
    }
    return sk;
}

// -------------------------- Ciphertext types --------------------------
struct CipherScalar {
    RnsPoly a;
    RnsPoly b;
    uint64_t correction = 1; // tracks accumulated q_last^{-1} mod p from mod-switches
};

struct CipherVector {
    RnsPoly a;
    vector<RnsPoly> b; // ell_prime polys, all sharing a
    uint64_t correction = 1;
};

struct TensorCt {
    RnsPoly d0;             // a1 * a2
    vector<RnsPoly> d1;     // a1 * b2[i]
    RnsPoly d2;             // a2 * b1
    vector<RnsPoly> d3;     // b1 * b2[i]
    uint64_t correction = 1;
};

// Result of mk-add: joint ciphertext under (s_1, s_2[0..ell'-1]).
//   Decrypt: out_msg[i] = b[i] - a1*s_1 - a2*s_2[i] = μ_1 + μ_2[i]  (mod p)
struct MkAddCt {
    RnsPoly a1;             // unchanged from ct_1
    RnsPoly a2;             // unchanged from ct_2 (shared across i)
    vector<RnsPoly> b;      // ell_prime polys; b[i] = b_1 + b_2[i]
    uint64_t correction = 1;
};

// -------------------------- Encryption --------------------------
static CipherScalar encrypt_scalar(const SecretKeyScalar& sk, const vector<uint64_t>& msg,
                                   const Params& params, Prng& prng) {
    if (msg.size() != params.N) throw std::runtime_error("msg size != N");
    const size_t L = params.q_chain.size();
    CipherScalar ct;
    ct.a = sample_uniform(params, L, prng);
    vector<int> e_coeffs = sample_noise_cbd(params.N, params.sigma_eta, prng);
    RnsPoly e_ntt = small_signed_to_rns_ntt(e_coeffs, params, L);
    RnsPoly as = mul(ct.a, sk.s_ntt, params);
    RnsPoly pe = mul_by_p(e_ntt, params);
    RnsPoly m_ntt = plain_to_rns_ntt(msg, params, L);
    ct.b = add(add(as, pe, params), m_ntt, params);
    ct.correction = 1;
    return ct;
}

static CipherVector encrypt_vector(const SecretKeyVector& sk,
                                   const vector<vector<uint64_t>>& msgs_m,
                                   const Params& params, Prng& prng) {
    if (msgs_m.size() != params.ell) throw std::runtime_error("msg vec len != ell");
    for (auto& m : msgs_m) if (m.size() != params.N) throw std::runtime_error("msg row size != N");

    const size_t L = params.q_chain.size();

    // Compute m' = T*m  in R_p^tau
    vector<vector<uint64_t>> mprime(params.tau, vector<uint64_t>(params.N, 0));
    for (size_t r = 0; r < params.tau; r++) {
        for (size_t c = 0; c < params.ell; c++) {
            auto prod = plain_mul(sk.T_coeffs[r][c], msgs_m[c], params);
            mprime[r] = plain_add(mprime[r], prod, params);
        }
    }
    // μ = m || m' of length ell_prime
    vector<vector<uint64_t>> mu(params.ell_prime);
    for (size_t i = 0; i < params.ell; i++) mu[i] = msgs_m[i];
    for (size_t j = 0; j < params.tau; j++) mu[params.ell + j] = mprime[j];

    CipherVector ct;
    ct.a = sample_uniform(params, L, prng);
    ct.b.reserve(params.ell_prime);
    for (size_t i = 0; i < params.ell_prime; i++) {
        vector<int> e_coeffs = sample_noise_cbd(params.N, params.sigma_eta, prng);
        RnsPoly e_ntt = small_signed_to_rns_ntt(e_coeffs, params, L);
        RnsPoly as = mul(ct.a, sk.s_ntt[i], params);
        RnsPoly pe = mul_by_p(e_ntt, params);
        RnsPoly mi = plain_to_rns_ntt(mu[i], params, L);
        ct.b.emplace_back(add(add(as, pe, params), mi, params));
    }
    ct.correction = 1;
    return ct;
}

// -------------------------- Multi-key Mul (one layer) --------------------------
// Output tensor: d0 = a1*a2; d1[i] = a1*b2[i]; d2 = a2*b1; d3[i] = b1*b2[i]
static TensorCt mkMul(const CipherScalar& ct1, const CipherVector& ct2, const Params& params) {
    if (ct1.a.num_primes != ct2.a.num_primes) throw std::runtime_error("level mismatch in mkMul");
    TensorCt out;
    out.d0 = mul(ct1.a, ct2.a, params);
    out.d2 = mul(ct2.a, ct1.b, params);
    const size_t ellp = ct2.b.size();
    out.d1.reserve(ellp);
    out.d3.reserve(ellp);
    for (size_t i = 0; i < ellp; i++) {
        out.d1.emplace_back(mul(ct1.a, ct2.b[i], params));
        out.d3.emplace_back(mul(ct1.b, ct2.b[i], params));
    }
    out.correction = multiply_uint_mod(ct1.correction, ct2.correction, params.p);
    return out;
}

// -------------------------- Multi-key Add --------------------------
// Adds a scalar ct (under sk_1) and a vector ct (under sk_2) componentwise. The shared "a"
// poly from each side is kept; the b polys are summed component-wise. Result is a joint
// ciphertext whose i-th component decrypts to μ_1 + μ_2[i] under (sk_1, sk_2[i]).
//
// Note: this does NOT preserve the sparsification check, since adding the same μ_1 to all
// components of μ_2 breaks the linear relation μ[ell:] == T·μ[:ell] in general. That is fine
// for benchmarking; if needed for a verified output, the caller can re-encrypt μ_1 as a vector.
static MkAddCt mkAdd(const CipherScalar& ct1, const CipherVector& ct2, const Params& params) {
    if (ct1.a.num_primes != ct2.a.num_primes) throw std::runtime_error("level mismatch in mkAdd");
    if (ct1.correction != ct2.correction) throw std::runtime_error("correction mismatch in mkAdd");
    MkAddCt out;
    out.a1 = ct1.a;
    out.a2 = ct2.a;
    out.b.reserve(ct2.b.size());
    for (size_t i = 0; i < ct2.b.size(); i++) {
        out.b.emplace_back(add(ct1.b, ct2.b[i], params));
    }
    out.correction = ct1.correction;
    return out;
}

// -------------------------- Mod switch wrappers --------------------------
static uint64_t inv_qlast_mod_p_at(const Params& params, size_t L) {
    uint64_t q_last = params.q_chain[L - 1].value();
    uint64_t inv;
    if (!try_invert_uint_mod(q_last, params.p, inv)) throw std::runtime_error("inv mod p failed");
    return inv;
}

static void mod_switch(CipherScalar& ct, const Params& params) {
    const size_t L = ct.a.num_primes;
    if (L <= 1) return;
    uint64_t corr = inv_qlast_mod_p_at(params, L);
    mod_switch_drop_last_poly(ct.a, params);
    mod_switch_drop_last_poly(ct.b, params);
    ct.correction = multiply_uint_mod(ct.correction, corr, params.p);
}

static void mod_switch(CipherVector& ct, const Params& params) {
    const size_t L = ct.a.num_primes;
    if (L <= 1) return;
    uint64_t corr = inv_qlast_mod_p_at(params, L);
    mod_switch_drop_last_poly(ct.a, params);
    for (auto& b : ct.b) mod_switch_drop_last_poly(b, params);
    ct.correction = multiply_uint_mod(ct.correction, corr, params.p);
}

static void mod_switch(TensorCt& ct, const Params& params) {
    const size_t L = ct.d0.num_primes;
    if (L <= 1) return;
    uint64_t corr = inv_qlast_mod_p_at(params, L);
    mod_switch_drop_last_poly(ct.d0, params);
    mod_switch_drop_last_poly(ct.d2, params);
    for (auto& d : ct.d1) mod_switch_drop_last_poly(d, params);
    for (auto& d : ct.d3) mod_switch_drop_last_poly(d, params);
    ct.correction = multiply_uint_mod(ct.correction, corr, params.p);
}

// -------------------------- Reduce RNS-NTT poly -> plaintext coeffs in [0,p) --------------------------
// Supports L=1 (single prime, 44-bit q') and L=2 (uint128 reconstruction) for the prototype.
static vector<uint64_t> reduce_to_plain(const RnsPoly& c, const Params& params, uint64_t correction) {
    const size_t L = c.num_primes;
    vector<uint64_t> out(params.N);

    if (L == 1) {
        vector<uint64_t> coeffs(params.N);
        std::memcpy(coeffs.data(), c.data.data(), params.N * sizeof(uint64_t));
        inverse_ntt_negacyclic_harvey(coeffs.data(), *params.q_ntt[0]);
        const uint64_t q = params.q_chain[0].value();
        const uint64_t q_half = q >> 1;
        const uint64_t p_val = params.p.value();
        const uint64_t q_mod_p = barrett_reduce_64(q, params.p);
        for (size_t i = 0; i < params.N; i++) {
            uint64_t v = coeffs[i];
            uint64_t v_mod_p = barrett_reduce_64(v, params.p);
            if (v > q_half) {
                // negative branch: signed value is v - q; mod p = v_mod_p - q_mod_p (signed)
                if (v_mod_p >= q_mod_p) out[i] = v_mod_p - q_mod_p;
                else out[i] = p_val - (q_mod_p - v_mod_p);
            } else {
                out[i] = v_mod_p;
            }
        }
    } else if (L == 2) {
        vector<uint64_t> c0(params.N), c1(params.N);
        std::memcpy(c0.data(), row(c, 0), params.N * sizeof(uint64_t));
        std::memcpy(c1.data(), row(c, 1), params.N * sizeof(uint64_t));
        inverse_ntt_negacyclic_harvey(c0.data(), *params.q_ntt[0]);
        inverse_ntt_negacyclic_harvey(c1.data(), *params.q_ntt[1]);
        const uint64_t q0 = params.q_chain[0].value();
        const uint64_t q1 = params.q_chain[1].value();
        const __uint128_t Q = __uint128_t(q0) * __uint128_t(q1);
        const __uint128_t Q_half = Q >> 1;
        uint64_t inv_q0_mod_q1;
        if (!try_invert_uint_mod(q0, params.q_chain[1], inv_q0_mod_q1))
            throw std::runtime_error("q0 not invertible mod q1");
        const int64_t p_val = (int64_t)params.p.value();
        for (size_t i = 0; i < params.N; i++) {
            const uint64_t r0 = c0[i], r1 = c1[i];
            uint64_t r0_mod_q1 = barrett_reduce_64(r0, params.q_chain[1]);
            uint64_t y = (r1 >= r0_mod_q1) ? (r1 - r0_mod_q1) : (q1 - (r0_mod_q1 - r1));
            y = multiply_uint_mod(y, inv_q0_mod_q1, params.q_chain[1]);
            __uint128_t x = __uint128_t(r0) + __uint128_t(q0) * __uint128_t(y);
            __int128_t xs = (x > Q_half) ? (__int128_t)x - (__int128_t)Q : (__int128_t)x;
            int64_t r = (int64_t)(xs % p_val);
            if (r < 0) r += p_val;
            out[i] = (uint64_t)r;
        }
    } else {
        throw std::runtime_error("reduce_to_plain: only L in {1,2} supported in this prototype; mod-switch first");
    }

    // Apply correction inverse: the decoded value equals (true_msg * correction) mod p
    if (correction != 1) {
        uint64_t inv_corr;
        if (!try_invert_uint_mod(correction, params.p, inv_corr))
            throw std::runtime_error("correction not invertible mod p");
        multiply_poly_scalar_coeffmod(out.data(), params.N, inv_corr, params.p, out.data());
    }
    return out;
}

// -------------------------- Decryption --------------------------
static vector<uint64_t> decrypt_scalar(const SecretKeyScalar& sk, const CipherScalar& ct, const Params& params) {
    const size_t L = ct.a.num_primes;
    RnsPoly s_lvl = truncate_to_level(sk.s_ntt, params, L);
    RnsPoly as = mul(ct.a, s_lvl, params);
    RnsPoly bma = sub(ct.b, as, params);
    return reduce_to_plain(bma, params, ct.correction);
}

// Add two CipherScalars with matching corrections. Plaintext-side result is μ_x + μ_y mod p.
static CipherScalar add_scalar(const CipherScalar& x, const CipherScalar& y, const Params& params) {
    if (x.a.num_primes != y.a.num_primes) throw std::runtime_error("level mismatch in add_scalar");
    if (x.correction != y.correction) throw std::runtime_error("correction mismatch in add_scalar");
    CipherScalar out;
    out.a = add(x.a, y.a, params);
    out.b = add(x.b, y.b, params);
    out.correction = x.correction;
    return out;
}

// Compute infinity-norm of the noise polynomial in (b - a*s) - μ_encoded, in coefficient form.
// Returns log2(max |e_i|) as a double; returns -1.0 if the decoded residue is inconsistent
// with the given μ (meaning decryption would fail and we can't pin down the noise).
// μ_in is the *expected* plaintext (after applying correction and reduction). Internally we
// compare against μ_encoded = (μ_in * ct.correction) mod p, recover e_signed s.t.
// (a*s + b)_centered = p*e_signed + μ_encoded, then return log2(max|e|).
static double noise_log2_scalar(const SecretKeyScalar& sk, const CipherScalar& ct,
                                const vector<uint64_t>& mu_in, const Params& params) {
    const size_t L = ct.a.num_primes;
    RnsPoly s_lvl = truncate_to_level(sk.s_ntt, params, L);
    RnsPoly as = mul(ct.a, s_lvl, params);
    RnsPoly bma = sub(ct.b, as, params);
    const uint64_t p_val = params.p.value();
    __int128_t max_abs = 0;
    auto pick_e = [&](__int128_t v) -> __int128_t {
        // v = p*e + mu_encoded, mu_encoded in [0, p). Try (v - mu_encoded)/p, else (v - (mu_encoded - p))/p
        return v;
    };
    if (L == 1) {
        vector<uint64_t> c(params.N);
        std::memcpy(c.data(), row(bma, 0), params.N * sizeof(uint64_t));
        inverse_ntt_negacyclic_harvey(c.data(), *params.q_ntt[0]);
        const uint64_t q = params.q_chain[0].value();
        const uint64_t q_half = q >> 1;
        for (size_t i = 0; i < params.N; i++) {
            __int128_t v = (c[i] > q_half) ? (__int128_t)c[i] - (__int128_t)q : (__int128_t)c[i];
            __int128_t mu_e = (__int128_t)((mu_in[i] * ct.correction) % p_val);
            __int128_t diff = v - mu_e;
            if (diff % (int64_t)p_val != 0) {
                __int128_t diff2 = v - (mu_e - (int64_t)p_val);
                if (diff2 % (int64_t)p_val != 0) return -1.0;
                diff = diff2;
            }
            __int128_t e = diff / (int64_t)p_val;
            if (e < 0) e = -e;
            if (e > max_abs) max_abs = e;
        }
    } else if (L == 2) {
        vector<uint64_t> c0(params.N), c1(params.N);
        std::memcpy(c0.data(), row(bma, 0), params.N * sizeof(uint64_t));
        std::memcpy(c1.data(), row(bma, 1), params.N * sizeof(uint64_t));
        inverse_ntt_negacyclic_harvey(c0.data(), *params.q_ntt[0]);
        inverse_ntt_negacyclic_harvey(c1.data(), *params.q_ntt[1]);
        const uint64_t q0 = params.q_chain[0].value();
        const uint64_t q1 = params.q_chain[1].value();
        const __uint128_t Q = __uint128_t(q0) * __uint128_t(q1);
        const __uint128_t Q_half = Q >> 1;
        uint64_t inv_q0_mod_q1;
        if (!try_invert_uint_mod(q0, params.q_chain[1], inv_q0_mod_q1))
            throw std::runtime_error("q0 not invertible mod q1");
        for (size_t i = 0; i < params.N; i++) {
            uint64_t r0_mod_q1 = barrett_reduce_64(c0[i], params.q_chain[1]);
            uint64_t y = (c1[i] >= r0_mod_q1) ? (c1[i] - r0_mod_q1) : (q1 - (r0_mod_q1 - c1[i]));
            y = multiply_uint_mod(y, inv_q0_mod_q1, params.q_chain[1]);
            __uint128_t x = __uint128_t(c0[i]) + __uint128_t(q0) * __uint128_t(y);
            __int128_t v = (x > Q_half) ? (__int128_t)x - (__int128_t)Q : (__int128_t)x;
            __int128_t mu_e = (__int128_t)((mu_in[i] * ct.correction) % p_val);
            __int128_t diff = v - mu_e;
            if (diff % (int64_t)p_val != 0) {
                __int128_t diff2 = v - (mu_e - (int64_t)p_val);
                if (diff2 % (int64_t)p_val != 0) return -1.0;
                diff = diff2;
            }
            __int128_t e = diff / (int64_t)p_val;
            if (e < 0) e = -e;
            if (e > max_abs) max_abs = e;
        }
    } else {
        throw std::runtime_error("noise_log2_scalar supports L in {1,2} only");
    }
    if (max_abs == 0) return 0.0;
    // log2 of a __int128
    int hi = 0;
    __uint128_t u = (__uint128_t)max_abs;
    while (u > 1) { u >>= 1; hi++; }
    // refine with bottom bits
    double frac = std::log2((double)(__uint128_t)max_abs / ((__uint128_t)1 << hi));
    return hi + frac;
    (void)pick_e;
}

static vector<vector<uint64_t>> decrypt_vector(const SecretKeyVector& sk, const CipherVector& ct,
                                                const Params& params, bool& check_passed) {
    const size_t L = ct.a.num_primes;
    vector<vector<uint64_t>> mu_full(params.ell_prime);
    for (size_t i = 0; i < params.ell_prime; i++) {
        RnsPoly s_lvl = truncate_to_level(sk.s_ntt[i], params, L);
        RnsPoly as = mul(ct.a, s_lvl, params);
        RnsPoly bma = sub(ct.b[i], as, params);
        mu_full[i] = reduce_to_plain(bma, params, ct.correction);
    }
    // Sparsification check: μ[ell+r] ?= sum_c T[r][c] · μ[c]
    check_passed = true;
    for (size_t r = 0; r < params.tau; r++) {
        vector<uint64_t> expected(params.N, 0);
        for (size_t c = 0; c < params.ell; c++) {
            auto prod = plain_mul(sk.T_coeffs[r][c], mu_full[c], params);
            expected = plain_add(expected, prod, params);
        }
        if (expected != mu_full[params.ell + r]) check_passed = false;
    }
    vector<vector<uint64_t>> out(params.ell);
    for (size_t i = 0; i < params.ell; i++) out[i] = mu_full[i];
    return out;
}

// μ[i] = d0·s1·s2[i] - d1[i]·s1 - d2·s2[i] + d3[i]   mod R_p
static vector<vector<uint64_t>> mkDec(const SecretKeyScalar& sk1, const SecretKeyVector& sk2,
                                      const TensorCt& tct, const Params& params, bool& check_passed) {
    const size_t L = tct.d0.num_primes;
    RnsPoly s1_lvl = truncate_to_level(sk1.s_ntt, params, L);

    vector<vector<uint64_t>> mu_full(params.ell_prime);
    for (size_t i = 0; i < params.ell_prime; i++) {
        RnsPoly s2i_lvl = truncate_to_level(sk2.s_ntt[i], params, L);

        RnsPoly t1 = mul(mul(tct.d0, s1_lvl, params), s2i_lvl, params); // d0*s1*s2[i]
        RnsPoly t2 = mul(tct.d1[i], s1_lvl, params);                    // d1[i]*s1
        RnsPoly t3 = mul(tct.d2, s2i_lvl, params);                      // d2*s2[i]
        RnsPoly r  = add(sub(sub(t1, t2, params), t3, params), tct.d3[i], params);
        mu_full[i] = reduce_to_plain(r, params, tct.correction);
    }
    check_passed = true;
    for (size_t r = 0; r < params.tau; r++) {
        vector<uint64_t> expected(params.N, 0);
        for (size_t c = 0; c < params.ell; c++) {
            auto prod = plain_mul(sk2.T_coeffs[r][c], mu_full[c], params);
            expected = plain_add(expected, prod, params);
        }
        if (expected != mu_full[params.ell + r]) check_passed = false;
    }
    vector<vector<uint64_t>> out(params.ell);
    for (size_t i = 0; i < params.ell; i++) out[i] = mu_full[i];
    return out;
}

// -------------------------- Helpers for reporting --------------------------
static int bit_length(uint64_t v) { int b = 0; while (v) { b++; v >>= 1; } return b; }
static size_t count_nonzero(const vector<int>& v) { size_t c = 0; for (int x : v) if (x) c++; return c; }

// -------------------------- Benchmark --------------------------
// Times mkMul and mkAdd at levels L=3 (full, ~130-bit q) and L=2 (after one mod-switch).
// Each operation is run `iters` times on the same input ciphertexts; we take the average.
static int run_benchmark(int iters = 100) {
    using namespace std::chrono;
    using std::cout;

    Params params;
    const uint64_t plain_modulus = (uint64_t(7) << 20) + 1;
    vector<int> q_bits = {44, 43, 43};
    params_init(params, plain_modulus, q_bits, /*N*/ 8192, /*ell*/ 4, /*tau*/ 1, /*hw*/ 32);

    int total_bits = 0;
    for (auto& q : params.q_chain) total_bits += bit_length(q.value());
    int last_bits = bit_length(params.q_chain.back().value());
    int bits_L2 = total_bits - last_bits;
    int bits_L1 = bit_length(params.q_chain.front().value());

    cout << "=== Benchmark: BGV multi-key ops, averaged over " << iters << " iterations ===\n";
    cout << "  N = " << params.N << ", p = " << plain_modulus
         << ",  ell = " << params.ell << ", tau = " << params.tau
         << ", ell' = " << params.ell_prime << ", hw = " << params.hw << "\n";
    cout << "  log q  @ L=3 = " << total_bits << " bits\n";
    cout << "  log q  @ L=2 = " << bits_L2    << " bits  (one mod-switch dropped " << last_bits << "-bit prime)\n";
    cout << "  log q' @ L=1 = " << bits_L1    << " bits  (for reference)\n\n";

    Prng prng(12345);

    SecretKeyScalar sk1 = keygen_scalar(params, prng);
    SecretKeyVector sk2 = keygen_vector(params, prng);

    vector<uint64_t> mu1(params.N);
    vector<vector<uint64_t>> m2(params.ell, vector<uint64_t>(params.N));
    for (size_t i = 0; i < params.N; i++) {
        mu1[i] = prng.uniform_below(plain_modulus);
        for (size_t c = 0; c < params.ell; c++) m2[c][i] = prng.uniform_below(plain_modulus);
    }

    // Fresh ciphertexts at L=3
    auto ct1_L3 = encrypt_scalar(sk1, mu1, params, prng);
    auto ct2_L3 = encrypt_vector(sk2, m2, params, prng);

    // L=2 copies via one mod-switch
    auto ct1_L2 = ct1_L3;  mod_switch(ct1_L2, params);
    auto ct2_L2 = ct2_L3;  mod_switch(ct2_L2, params);

    auto bench_mul = [&](const CipherScalar& a, const CipherVector& b, const std::string& label) {
        // Warm-up (1 untimed call so allocator/page-cache settles)
        { auto w = mkMul(a, b, params); (void)w.d0.data[0]; }
        uint64_t sink = 0;
        auto t0 = steady_clock::now();
        for (int i = 0; i < iters; i++) {
            auto r = mkMul(a, b, params);
            sink ^= r.d0.data[0] ^ r.d3.back().data[0];  // touch result to avoid DCE
        }
        auto t1 = steady_clock::now();
        double avg_ms = duration_cast<nanoseconds>(t1 - t0).count() / double(iters) / 1.0e6;
        cout << "  " << std::left << std::setw(14) << label
             << "  avg " << std::fixed << std::setprecision(3) << std::setw(8) << avg_ms << " ms"
             << "  (over " << iters << " iters)\n";
        (void)sink;
    };

    auto bench_add = [&](const CipherScalar& a, const CipherVector& b, const std::string& label) {
        { auto w = mkAdd(a, b, params); (void)w.b[0].data[0]; }
        uint64_t sink = 0;
        auto t0 = steady_clock::now();
        for (int i = 0; i < iters; i++) {
            auto r = mkAdd(a, b, params);
            sink ^= r.b.front().data[0] ^ r.b.back().data[0];
        }
        auto t1 = steady_clock::now();
        double avg_ms = duration_cast<nanoseconds>(t1 - t0).count() / double(iters) / 1.0e6;
        cout << "  " << std::left << std::setw(14) << label
             << "  avg " << std::fixed << std::setprecision(3) << std::setw(8) << avg_ms << " ms"
             << "  (over " << iters << " iters)\n";
        (void)sink;
    };

    bench_mul(ct1_L3, ct2_L3, "mkMul @ L=3");
    bench_add(ct1_L3, ct2_L3, "mkAdd @ L=3");
    bench_mul(ct1_L2, ct2_L2, "mkMul @ L=2");
    bench_add(ct1_L2, ct2_L2, "mkAdd @ L=2");

    return 0;
}

// Synthetic ciphertext statistically equivalent to summing K independent fresh encryptions
// of μ. Doing 2^21 real encryptions is too slow, but their sum is just (Σa_i, Σa_i·s + p·Σe_i + Kμ):
//   - Σa_i (sum of K independent uniforms mod q) is itself uniform mod q.
//   - Σe_i has each coefficient ~ Normal(0, K·η/2) by CLT (η = CBD parameter).
// So we sample one uniform a and one Gaussian-noise e directly.
static CipherScalar simulate_sum_K_fresh(const SecretKeyScalar& sk, const vector<uint64_t>& mu,
                                         uint64_t K, const Params& params, Prng& prng) {
    const size_t L = params.q_chain.size();
    CipherScalar ct;
    ct.a = sample_uniform(params, L, prng);

    const double sigma_sum = std::sqrt((double)K * (double)params.sigma_eta / 2.0);
    std::mt19937_64 rng(prng.u64());
    std::normal_distribution<double> nd(0.0, sigma_sum);
    vector<int64_t> e_signed(params.N);
    for (size_t i = 0; i < params.N; i++) {
        e_signed[i] = (int64_t)std::llround(nd(rng));
    }

    // RNS-NTT-encode e_signed (per prime).
    RnsPoly e_ntt = make_zero(params, L);
    for (size_t k = 0; k < L; k++) {
        uint64_t q = params.q_chain[k].value();
        uint64_t* dst = row(e_ntt, k);
        for (size_t i = 0; i < params.N; i++) {
            int64_t v = e_signed[i];
            if (v < 0) {
                int64_t r = (-v) % (int64_t)q;
                dst[i] = (r == 0) ? 0 : (q - (uint64_t)r);
            } else {
                dst[i] = (uint64_t)(v % (int64_t)q);
            }
        }
        ntt_negacyclic_harvey(dst, *params.q_ntt[k]);
    }

    // Σμ_i = K·μ. We can compute K mod p first.
    vector<uint64_t> Kmu(params.N);
    uint64_t K_mod_p = K % params.p.value();
    for (size_t i = 0; i < params.N; i++) {
        Kmu[i] = multiply_uint_mod(mu[i], K_mod_p, params.p);
    }

    RnsPoly as = mul(ct.a, sk.s_ntt, params);
    RnsPoly pe = mul_by_p(e_ntt, params);
    RnsPoly mu_ntt = plain_to_rns_ntt(Kmu, params, L);
    ct.b = add(add(as, pe, params), mu_ntt, params);
    ct.correction = 1;
    return ct;
}

// Add p · uniform([-2^smudge_log2, 2^smudge_log2]) to `target` (an R_q poly in NTT form).
// The smudge enters the noise side (multiplied by p), so it doesn't change the plaintext.
static void add_smudge_poly(RnsPoly& target, int smudge_log2, const Params& params, Prng& prng) {
    if (smudge_log2 < 1 || smudge_log2 > 62) throw std::runtime_error("smudge_log2 out of range");
    const size_t L = target.num_primes;
    const int64_t bound = (int64_t)1 << smudge_log2;
    std::mt19937_64 rng(prng.u64());
    std::uniform_int_distribution<int64_t> uid(-(bound - 1), bound - 1);

    vector<int64_t> s_signed(params.N);
    for (size_t i = 0; i < params.N; i++) s_signed[i] = uid(rng);

    RnsPoly sm = make_zero(params, L);
    for (size_t k = 0; k < L; k++) {
        uint64_t q = params.q_chain[k].value();
        uint64_t* dst = row(sm, k);
        for (size_t i = 0; i < params.N; i++) {
            int64_t v = s_signed[i];
            if (v < 0) {
                int64_t r = (-v) % (int64_t)q;
                dst[i] = (r == 0) ? 0 : (q - (uint64_t)r);
            } else {
                dst[i] = (uint64_t)(v % (int64_t)q);
            }
        }
        ntt_negacyclic_harvey(dst, *params.q_ntt[k]);
    }
    RnsPoly p_sm = mul_by_p(sm, params);
    target = add(target, p_sm, params);
}

// Add a uniform smudging noise of magnitude 2^smudge_log2 (coefficient-wise, signed) to ct's b
// polynomial. The smudge enters on the noise side (multiplied by p), so it adds to the BGV noise
// term and leaves the encoded message unchanged. After this call, ||e||_inf ≈ 2^smudge_log2 if
// smudge dominates the pre-existing noise.
static void add_smudge(CipherScalar& ct, int smudge_log2, const Params& params, Prng& prng) {
    if (smudge_log2 < 1 || smudge_log2 > 62) throw std::runtime_error("smudge_log2 out of range");
    const size_t L = ct.a.num_primes;
    const int64_t bound = (int64_t)1 << smudge_log2;
    std::mt19937_64 rng(prng.u64());
    std::uniform_int_distribution<int64_t> uid(-(bound - 1), bound - 1);

    RnsPoly sm_ntt = make_zero(params, L);
    for (size_t k = 0; k < L; k++) {
        uint64_t q = params.q_chain[k].value();
        uint64_t* dst = row(sm_ntt, k);
        // Sample fresh per coord (same per-prime by drawing one signed value and reducing mod q).
    }
    // Sample once, then reduce mod each q.
    vector<int64_t> s_signed(params.N);
    for (size_t i = 0; i < params.N; i++) s_signed[i] = uid(rng);

    for (size_t k = 0; k < L; k++) {
        uint64_t q = params.q_chain[k].value();
        uint64_t* dst = row(sm_ntt, k);
        for (size_t i = 0; i < params.N; i++) {
            int64_t v = s_signed[i];
            if (v < 0) {
                int64_t r = (-v) % (int64_t)q;
                dst[i] = (r == 0) ? 0 : (q - (uint64_t)r);
            } else {
                dst[i] = (uint64_t)(v % (int64_t)q);
            }
        }
        ntt_negacyclic_harvey(dst, *params.q_ntt[k]);
    }
    RnsPoly p_sm = mul_by_p(sm_ntt, params);
    ct.b = add(ct.b, p_sm, params);
}

// -------------------------- Noise test: 2^21 additions at L=2 --------------------------
// Strategy: at L=2 we have ct_0 (fresh-mod-switched). We then double 21 times:
//   ct_{k+1} = ct_k + ct_k       (one add operation per level)
// After 21 levels, the resulting ct represents 2^{21} · μ mod p — equivalent to a binary tree
// of 2^{21} leaves, all equal to ct_0. Total ct-additions performed: 21. The noise grows
// deterministically by a factor of 2 each level, so ||e_final||_∞ = 2^{21} · ||e_0||_∞ exactly.
//
// (If we used 2^{21} *fresh-independent* ciphertexts at the leaves, the noise std-dev would
//  grow as sqrt(2^{21}) · σ_fresh ≈ 1448 · σ, much smaller. We report both for comparison.)
static int run_add_noise_test(int q_last_bits = 44) {
    using std::cout;
    using std::endl;

    Params params;
    const uint64_t plain_modulus = (uint64_t(7) << 20) + 1;
    vector<int> q_bits = {q_last_bits, 43, 43};
    params_init(params, plain_modulus, q_bits, /*N*/ 8192, /*ell*/ 4, /*tau*/ 1, /*hw*/ 32);

    int total_bits = 0;
    int last_bits = bit_length(params.q_chain.back().value());
    for (auto& q : params.q_chain) total_bits += bit_length(q.value());
    int bits_L2 = total_bits - last_bits;
    int bits_L1 = bit_length(params.q_chain.front().value());

    cout << "=== Noise test: 2^21 additions at L=2 (binary tree of doublings) ===\n";
    cout << "  N = " << params.N << ", p = " << plain_modulus
         << " (log2 p ≈ " << std::log2((double)plain_modulus) << " bits)\n";
    cout << "  q chain (drop order right-to-left):";
    for (auto& q : params.q_chain) cout << "  " << bit_length(q.value()) << "b";
    cout << "\n  log q @ L=3 = " << total_bits << ",  L=2 = " << bits_L2 << ",  L=1 = " << bits_L1 << "\n\n";

    Prng prng(424242);

    SecretKeyScalar sk1 = keygen_scalar(params, prng);

    vector<uint64_t> mu(params.N);
    for (size_t i = 0; i < params.N; i++) mu[i] = prng.uniform_below(plain_modulus);

    // Encrypt at L=3, then mod-switch once to L=2. (Noise measurement runs at L=2; the L=3
    // reconstruction would need 130-bit CRT which overflows __int128 by a hair.)
    auto ct = encrypt_scalar(sk1, mu, params, prng);
    mod_switch(ct, params);  // L=3 → L=2
    double n0 = noise_log2_scalar(sk1, ct, mu, params);
    cout << std::fixed << std::setprecision(2)
         << "[L=2 after mod-switch] noise log2(||e||_inf) = " << n0
         << "  (budget ≈ " << bits_L2 - std::log2((double)plain_modulus) - 1.0 << " bits)\n\n";

    // Verify ct still decrypts correctly before adds.
    {
        auto dec = decrypt_scalar(sk1, ct, params);
        cout << "[pre-add decrypt] " << (dec == mu ? "OK" : "FAIL") << "\n";
    }

    // 21 doublings — binary tree of 2^21 leaves, all = ct.
    auto t0 = std::chrono::steady_clock::now();
    for (int k = 1; k <= 21; k++) {
        ct = add_scalar(ct, ct, params);
    }
    auto t1 = std::chrono::steady_clock::now();
    cout << "\n[after 21 doublings (2^21 effective adds)] "
         << std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() << " us total\n";

    // Expected plaintext: (2^21 · μ) mod p
    vector<uint64_t> expected(params.N);
    const uint64_t scale = (uint64_t(1) << 21) % plain_modulus;
    for (size_t i = 0; i < params.N; i++) {
        expected[i] = multiply_uint_mod(mu[i], scale, params.p);
    }

    double n_final = noise_log2_scalar(sk1, ct, expected, params);
    cout << "  noise log2(||e||_inf) = " << n_final;
    if (n_final >= 0) {
        double budget = bits_L2 - std::log2((double)plain_modulus) - 1.0;
        cout << "   (budget at L=2 ≈ " << budget << " bits, remaining ≈ "
             << (budget - n_final) << " bits)\n";
    } else {
        cout << "  (CT inconsistent — decryption would fail before reaching here)\n";
    }

    // Decrypt and check correctness.
    auto dec = decrypt_scalar(sk1, ct, params);
    bool ok = (dec == expected);
    cout << "\n[decrypt] " << (ok ? "OK" : "FAIL") << " — result " << (ok ? "matches" : "DOES NOT MATCH")
         << " (2^21 · μ) mod p\n";

    // Report theoretical noise for the "fresh-independent leaves" alternative.
    cout << "\n[theoretical noise comparison at L=2 after 2^21 adds]\n"
         << "  doubling (worst-case, what we just did):   ||e||_inf ≈ 2^21 · ||e_0||_inf"
         << "  → log2 ≈ " << (n0 + 21.0) << " bits\n"
         << "  fresh-independent leaves (random walk):    ||e||_inf ≈ sqrt(2^21) · σ_fresh"
         << "  → log2 ≈ " << (n0 + 10.5) << " bits (approx)\n";

    return ok ? 0 : 1;
}

// -------------------------- Smudge + mod-switch pipeline test --------------------------
// Pipeline:
//   1) Synthesize a ct at L=3 statistically equivalent to summing K=2^21 fresh encryptions of μ.
//   2) Mod-switch L=3 → L=2 (drops q2).
//   3) Add a uniform smudge of magnitude 2^{smudge_log2} (uses ~smudge_log2 - log2|e_pre| budget).
//   4) Mod-switch L=2 → L=1 (drops q1, dividing the noise by ~q1).
//   5) Decrypt at L=1 (single prime q0 of `q0_bits` bits). Check result == (K · μ) mod p.
// Also report final ciphertext size at L=1 = 2·N·q0_bits / 8 bytes.
static int run_smudge_test(int q0_bits = 40, int smudge_log2 = 50, int K_log2 = 21) {
    using std::cout;

    Params params;
    const uint64_t plain_modulus = (uint64_t(7) << 20) + 1;
    vector<int> q_bits = {q0_bits, 43, 43};
    params_init(params, plain_modulus, q_bits, /*N*/ 8192, /*ell*/ 4, /*tau*/ 1, /*hw*/ 32);

    int last_bits = bit_length(params.q_chain.back().value());
    int bits_L3 = 0;
    for (auto& q : params.q_chain) bits_L3 += bit_length(q.value());
    int bits_L2 = bits_L3 - last_bits;
    int q1_bits = bit_length(params.q_chain[1].value());
    int q0_actual = bit_length(params.q_chain[0].value());

    cout << "=== Smudge pipeline test ===\n"
         << "  N = " << params.N << ", p = " << plain_modulus
         << " (log2 p ≈ " << std::log2((double)plain_modulus) << " b)\n"
         << "  q chain (drop right-to-left):  q0=" << q0_actual << "b  q1="
         << q1_bits << "b  q2=" << last_bits << "b   (total " << bits_L3 << "b)\n"
         << "  pipeline parameters: K = 2^" << K_log2 << ", smudge = 2^" << smudge_log2 << "\n\n";

    Prng prng(7777);
    SecretKeyScalar sk = keygen_scalar(params, prng);
    vector<uint64_t> mu(params.N);
    for (size_t i = 0; i < params.N; i++) mu[i] = prng.uniform_below(plain_modulus);

    // Step 1: synthesize sum of K = 2^K_log2 fresh ciphertexts at L=3.
    const uint64_t K = (uint64_t(1) << K_log2);
    auto ct = simulate_sum_K_fresh(sk, mu, K, params, prng);

    // Step 2: mod-switch L=3 → L=2.
    mod_switch(ct, params);
    vector<uint64_t> K_mu(params.N);
    uint64_t K_mod_p = K % plain_modulus;
    for (size_t i = 0; i < params.N; i++) K_mu[i] = multiply_uint_mod(mu[i], K_mod_p, params.p);

    double n_post_adds = noise_log2_scalar(sk, ct, K_mu, params);
    cout << std::fixed << std::setprecision(2)
         << "[L=2 after 2^" << K_log2 << " simulated independent adds]\n"
         << "    log2||e||_inf = " << n_post_adds
         << "   (budget at L=2 ≈ " << bits_L2 - std::log2((double)plain_modulus) - 1.0 << "b)\n";
    {
        auto dec = decrypt_scalar(sk, ct, params);
        cout << "    decrypt: " << (dec == K_mu ? "OK" : "FAIL") << "\n";
    }

    // Step 3: smudge.
    add_smudge(ct, smudge_log2, params, prng);
    double n_post_smudge = noise_log2_scalar(sk, ct, K_mu, params);
    cout << "\n[L=2 after smudging with 2^" << smudge_log2 << "]\n"
         << "    log2||e||_inf = " << n_post_smudge
         << "   (budget remaining ≈ " << (bits_L2 - std::log2((double)plain_modulus) - 1.0 - n_post_smudge) << "b)\n";
    {
        auto dec = decrypt_scalar(sk, ct, params);
        cout << "    decrypt: " << (dec == K_mu ? "OK" : "FAIL") << "\n";
    }

    // Step 4: mod-switch L=2 → L=1.
    mod_switch(ct, params);
    double n_post_modswitch = noise_log2_scalar(sk, ct, K_mu, params);
    cout << "\n[L=1 after mod-switch (dropped q1, " << q1_bits << "b)]\n"
         << "    log2||e||_inf = " << n_post_modswitch
         << "   (noise scaled by ~1/q1, so dropped by ~" << q1_bits << "b)\n"
         << "    L=1 budget ≈ " << q0_actual - std::log2((double)plain_modulus) - 1.0 << "b\n";

    // Step 5: decrypt and check.
    auto dec = decrypt_scalar(sk, ct, params);
    bool ok = (dec == K_mu);
    cout << "    decrypt: " << (ok ? "OK" : "FAIL")
         << " — result " << (ok ? "matches" : "DOES NOT MATCH") << " (2^" << K_log2 << " · μ) mod p\n";

    // Ciphertext size at L=1.
    // 2 polys (a, b) × N coefficients, each coefficient fits in q0_actual bits.
    double bytes = 2.0 * params.N * q0_actual / 8.0;
    cout << "\n[final ciphertext at L=1]"
         << "\n    components:     2 R_q polys (a, b)"
         << "\n    coefficients:   " << params.N << " per poly"
         << "\n    coeff width:    " << q0_actual << " bits"
         << "\n    total size:     " << (size_t)bytes << " bytes (" << bytes / 1024.0 << " KiB)\n";

    return ok ? 0 : 1;
}

// -------------------------- Multi-key smudge pipeline STRESS test --------------------------
// Repeats the multi-key smudge pipeline `iters` times with independent randomness (keys,
// messages, encryption noise, smudge), counts decryption failures at L=2 (after smudge,
// before final mod-switch) and at L=1 (after final mod-switch).
static int run_smudge_mk_stress(int q0_bits, int smudge_log2, int iters) {
    using std::cout;
    Params params;
    const uint64_t plain_modulus = (uint64_t(7) << 20) + 1;
    vector<int> q_bits = {q0_bits, 43, 43};
    params_init(params, plain_modulus, q_bits, /*N*/ 8192, /*ell*/ 4, /*tau*/ 1, /*hw*/ 32);

    int q0_actual = bit_length(params.q_chain[0].value());
    int q1_bits   = bit_length(params.q_chain[1].value());
    int q2_bits   = bit_length(params.q_chain[2].value());

    cout << "=== Multi-key smudge pipeline STRESS test ===\n"
         << "  N = " << params.N << ", p = " << plain_modulus
         << ", ell' = " << params.ell_prime << "\n"
         << "  q chain: q0=" << q0_actual << "b  q1=" << q1_bits << "b  q2=" << q2_bits << "b\n"
         << "  smudge magnitude = 2^" << smudge_log2 << " (on each d3[i])\n"
         << "  iterations:    " << iters << "\n\n";

    int ok_total = 0;
    int fail_L2 = 0, fail_L1_only = 0;
    int sparsi_fail_L1 = 0;
    int worst_progress_print = -1;
    auto t0 = std::chrono::steady_clock::now();

    for (int it = 0; it < iters; it++) {
        Prng prng((uint64_t)it * 1000003ULL + 7777ULL);

        SecretKeyScalar sk1 = keygen_scalar(params, prng);
        SecretKeyVector sk2 = keygen_vector(params, prng);

        vector<uint64_t> mu1(params.N);
        vector<vector<uint64_t>> m2(params.ell, vector<uint64_t>(params.N));
        for (size_t i = 0; i < params.N; i++) {
            mu1[i] = prng.uniform_below(plain_modulus);
            for (size_t c = 0; c < params.ell; c++) m2[c][i] = prng.uniform_below(plain_modulus);
        }

        auto ct1 = encrypt_scalar(sk1, mu1, params, prng);
        auto ct2 = encrypt_vector(sk2, m2, params, prng);

        auto tct = mkMul(ct1, ct2, params);
        mod_switch(tct, params);  // L=3 → L=2

        for (size_t i = 0; i < tct.d3.size(); i++) {
            add_smudge_poly(tct.d3[i], smudge_log2, params, prng);
        }

        bool ok_chk;
        auto dec_L2 = mkDec(sk1, sk2, tct, params, ok_chk);
        bool match_L2 = true;
        for (size_t i = 0; i < params.ell; i++) {
            auto expected = plain_mul(mu1, m2[i], params);
            if (dec_L2[i] != expected) { match_L2 = false; break; }
        }

        mod_switch(tct, params);  // L=2 → L=1
        auto dec_L1 = mkDec(sk1, sk2, tct, params, ok_chk);
        bool match_L1 = true;
        for (size_t i = 0; i < params.ell; i++) {
            auto expected = plain_mul(mu1, m2[i], params);
            if (dec_L1[i] != expected) { match_L1 = false; break; }
        }

        if (match_L1) ok_total++;
        else if (!match_L2) fail_L2++;
        else fail_L1_only++;

        if (match_L1 && !ok_chk) sparsi_fail_L1++;

        int pct = (it + 1) * 100 / iters;
        if (pct >= worst_progress_print + 10) {
            auto t1 = std::chrono::steady_clock::now();
            double sec = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() / 1000.0;
            cout << "  [progress] " << (it + 1) << "/" << iters
                 << " (" << pct << "%)  ok=" << ok_total
                 << "  fail_L2=" << fail_L2 << "  fail_L1_only=" << fail_L1_only
                 << "   elapsed " << std::fixed << std::setprecision(1) << sec << "s\n";
            worst_progress_print = pct;
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    double sec = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() / 1000.0;

    cout << "\n=== Summary ===\n"
         << "  iterations:                  " << iters << "\n"
         << "  ok (final L=1 decrypt):      " << ok_total
         << "  (" << std::fixed << std::setprecision(3) << 100.0 * ok_total / iters << "%)\n"
         << "  failed at L=2 (post-smudge): " << fail_L2 << "\n"
         << "  failed at L=1 only:          " << fail_L1_only << "\n"
         << "  sparsification check fails:  " << sparsi_fail_L1 << "\n"
         << "  total runtime:               " << sec << "s ("
         << (sec * 1000.0 / iters) << " ms/iter)\n";

    return (ok_total == iters) ? 0 : 1;
}

// -------------------------- Multi-key smudge pipeline test --------------------------
// Pipeline (true multi-key, not scalar):
//   1) Encrypt ct_1 (scalar, BGV) and ct_2 (vector, ell'=5) at L=3.
//   2) mkMul → TensorCt with 1 + ell' + 1 + ell' = 2 + 2·ell' = 12 R_q-poly components.
//   3) Mod-switch tensor L=3 → L=2.
//   4) Smudge d3[i] for each i ∈ [ell'] with magnitude 2^smudge_log2. d3 enters mkDec without
//      any s-factor (the decryption is d0·s1·s2[i] - d1[i]·s1 - d2·s2[i] + d3[i]), so smudge
//      lands directly in the noise term. Smudging d0/d1/d2 also works but gets amplified by
//      ||s_1·s_2[i]||, so we'd need a smaller smudge magnitude there for the same effect.
//   5) Mod-switch tensor L=2 → L=1 with a small final prime q0.
//   6) mkDec; verify each slot decrypts to (μ_1 · μ_2[i]) mod p.
//   7) Report total ct size = 12 · N · q0_bits / 8.
static int run_smudge_mk_test(int q0_bits, int smudge_log2) {
    using std::cout;
    Params params;
    const uint64_t plain_modulus = (uint64_t(7) << 20) + 1;
    vector<int> q_bits = {q0_bits, 43, 43};
    params_init(params, plain_modulus, q_bits, /*N*/ 8192, /*ell*/ 4, /*tau*/ 1, /*hw*/ 32);

    int q0_actual = bit_length(params.q_chain[0].value());
    int q1_bits   = bit_length(params.q_chain[1].value());
    int q2_bits   = bit_length(params.q_chain[2].value());
    int bits_L3   = q0_actual + q1_bits + q2_bits;

    cout << "=== Multi-key smudge pipeline test (true mkMul tensor) ===\n"
         << "  N = " << params.N << ", p = " << plain_modulus
         << ", ell = " << params.ell << ", tau = " << params.tau
         << ", ell' = " << params.ell_prime << "\n"
         << "  q chain (drop right→left):  q0=" << q0_actual << "b  q1="
         << q1_bits << "b  q2=" << q2_bits << "b  (total " << bits_L3 << "b)\n"
         << "  smudge magnitude = 2^" << smudge_log2 << " (on d3[i] only)\n\n";

    Prng prng(7777);
    SecretKeyScalar sk1 = keygen_scalar(params, prng);
    SecretKeyVector sk2 = keygen_vector(params, prng);

    vector<uint64_t> mu1(params.N);
    vector<vector<uint64_t>> m2(params.ell, vector<uint64_t>(params.N));
    for (size_t i = 0; i < params.N; i++) {
        mu1[i] = prng.uniform_below(plain_modulus);
        for (size_t c = 0; c < params.ell; c++) m2[c][i] = prng.uniform_below(plain_modulus);
    }

    auto ct1 = encrypt_scalar(sk1, mu1, params, prng);
    auto ct2 = encrypt_vector(sk2, m2, params, prng);

    // mkMul at L=3 -> tensor with 12 components.
    auto tct = mkMul(ct1, ct2, params);
    const size_t tensor_polys = 1 + tct.d1.size() + 1 + tct.d3.size();
    cout << "[mkMul]  tensor at L=3:  1 d0 + " << tct.d1.size()
         << " d1 + 1 d2 + " << tct.d3.size() << " d3 = " << tensor_polys << " R_q polys\n";

    // Step 3: mod-switch to L=2.
    mod_switch(tct, params);
    cout << "[mod-switch] L=3 → L=2  (dropped " << q2_bits << "-bit prime)\n";
    {
        bool ok;
        auto dec = mkDec(sk1, sk2, tct, params, ok);
        bool match = true;
        for (size_t i = 0; i < params.ell; i++) {
            auto expected = plain_mul(mu1, m2[i], params);
            if (dec[i] != expected) match = false;
        }
        cout << "    pre-smudge mkDec @ L=2:  " << (match ? "OK" : "FAIL")
             << ", sparsification check: " << (ok ? "passed" : "FAILED") << "\n";
    }

    // Step 4: smudge d3[i] for each i.
    for (size_t i = 0; i < tct.d3.size(); i++) {
        add_smudge_poly(tct.d3[i], smudge_log2, params, prng);
    }
    {
        bool ok;
        auto dec = mkDec(sk1, sk2, tct, params, ok);
        bool match = true;
        for (size_t i = 0; i < params.ell; i++) {
            auto expected = plain_mul(mu1, m2[i], params);
            if (dec[i] != expected) match = false;
        }
        cout << "[smudge]      added 2^" << smudge_log2 << " to each d3[i] (" << tct.d3.size() << " components)\n"
             << "    post-smudge mkDec @ L=2: " << (match ? "OK" : "FAIL") << "\n";
    }

    // Step 5: mod-switch to L=1.
    mod_switch(tct, params);
    cout << "[mod-switch] L=2 → L=1  (dropped " << q1_bits << "-bit prime)\n";

    // Step 6: mkDec at L=1 and verify.
    bool ok;
    auto dec = mkDec(sk1, sk2, tct, params, ok);
    bool match = true;
    for (size_t i = 0; i < params.ell; i++) {
        auto expected = plain_mul(mu1, m2[i], params);
        if (dec[i] != expected) match = false;
    }
    cout << "    final mkDec @ L=1:       " << (match ? "OK" : "FAIL")
         << ", sparsification: " << (ok ? "passed" : "FAILED")
         << "  (final prime q0 = " << q0_actual << "b)\n";

    // Step 7: report size.
    double bytes = (double)tensor_polys * params.N * q0_actual / 8.0;
    cout << "\n[final ciphertext at L=1]"
         << "\n    components:   " << tensor_polys << " R_q polys (= 2 + 2·ell')"
         << "\n    coefficients: " << params.N << " per poly"
         << "\n    coeff width:  " << q0_actual << " bits"
         << "\n    total size:   " << (size_t)bytes << " bytes ("
         << bytes / 1024.0 << " KiB)\n";

    return match ? 0 : 1;
}

} // namespace mk

// ============================== main ==============================
int main(int argc, char* argv[]) {
    using namespace mk;
    using std::cout;
    using std::endl;

    // CLI dispatch: `./prototype bench [iters]` runs the benchmark; otherwise run the demo.
    if (argc >= 2) {
        std::string cmd = argv[1];
        if (cmd == "bench" || cmd == "--bench" || cmd == "--benchmark") {
            int iters = 100;
            if (argc >= 3) iters = std::max(1, std::atoi(argv[2]));
            return run_benchmark(iters);
        }
        if (cmd == "noise-add") {
            int last_bits = (argc >= 3) ? std::atoi(argv[2]) : 44;
            return run_add_noise_test(last_bits);
        }
        if (cmd == "smudge") {
            int q0_bits = (argc >= 3) ? std::atoi(argv[2]) : 40;
            int smudge = (argc >= 4) ? std::atoi(argv[3]) : 50;
            int K_log2 = (argc >= 5) ? std::atoi(argv[4]) : 21;
            return run_smudge_test(q0_bits, smudge, K_log2);
        }
        if (cmd == "smudge-mk") {
            int q0_bits = (argc >= 3) ? std::atoi(argv[2]) : 40;
            int smudge = (argc >= 4) ? std::atoi(argv[3]) : 50;
            return run_smudge_mk_test(q0_bits, smudge);
        }
        if (cmd == "smudge-mk-stress") {
            int q0_bits = (argc >= 3) ? std::atoi(argv[2]) : 32;
            int smudge = (argc >= 4) ? std::atoi(argv[3]) : 50;
            int iters = (argc >= 5) ? std::atoi(argv[4]) : 1000;
            return run_smudge_mk_stress(q0_bits, smudge, iters);
        }
    }

    cout << "=== Multi-key RLWE Vector Encryption — BGV prototype (one mul layer) ===\n";

    Params params;
    const uint64_t plain_modulus = (uint64_t(7) << 20) + 1; // 7340033 = 7*2^20 + 1
    // Three primes ≈ {44, 43, 43} bits, total ≈ 130 bits; manual mod-switch twice -> 1 prime ≈ 44 bits.
    vector<int> q_bits = {44, 43, 43};
    params_init(params, plain_modulus, q_bits, /*N*/ 8192, /*ell*/ 4, /*tau*/ 1, /*hw*/ 32);

    cout << "  N      = " << params.N << "\n";
    cout << "  p      = " << plain_modulus << "  (= 7*2^20 + 1)\n";
    int total_bits = 0;
    cout << "  q chain (in mod-switch drop order from the right):\n";
    for (size_t i = 0; i < params.q_chain.size(); i++) {
        int bits = bit_length(params.q_chain[i].value());
        total_bits += bits;
        cout << "    q[" << i << "] = " << params.q_chain[i].value() << "   (" << bits << " bits)\n";
    }
    cout << "  total log q = " << total_bits << " bits\n";
    cout << "  ell = " << params.ell << ", tau = " << params.tau
         << ", ell' = " << params.ell_prime << ", hw = " << params.hw << "\n\n";

    Prng prng(12345);

    auto t0 = std::chrono::steady_clock::now();
    SecretKeyScalar sk1 = keygen_scalar(params, prng);
    SecretKeyVector sk2 = keygen_vector(params, prng);
    auto t1 = std::chrono::steady_clock::now();

    cout << "[keygen] HW(s1) = " << count_nonzero(sk1.s_coeffs) << " (expected " << params.hw << ")\n";
    for (size_t i = 0; i < sk2.s_coeffs.size(); i++) {
        cout << "[keygen] HW(s2[" << i << "]) = " << count_nonzero(sk2.s_coeffs[i])
             << " (expected " << params.hw << ")\n";
    }
    cout << "[keygen] " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms\n\n";

    // Test messages
    vector<uint64_t> mu1(params.N);
    vector<vector<uint64_t>> m2(params.ell, vector<uint64_t>(params.N));
    for (size_t i = 0; i < params.N; i++) {
        mu1[i] = prng.uniform_below(plain_modulus);
        for (size_t c = 0; c < params.ell; c++) m2[c][i] = prng.uniform_below(plain_modulus);
    }

    // ---------- Encrypt ----------
    auto t_enc0 = std::chrono::steady_clock::now();
    auto ct1 = encrypt_scalar(sk1, mu1, params, prng);
    auto ct2 = encrypt_vector(sk2, m2, params, prng);
    auto t_enc1 = std::chrono::steady_clock::now();
    cout << "[encrypt] both ciphertexts at L=" << ct1.a.num_primes << " primes, "
         << std::chrono::duration_cast<std::chrono::milliseconds>(t_enc1 - t_enc0).count() << " ms\n";

    // ---------- Sanity: decrypt freshly-encrypted ciphertexts after 1 mod-switch (so L=2) ----------
    {
        auto ct1s = ct1;
        auto ct2s = ct2;
        mod_switch(ct1s, params);
        mod_switch(ct2s, params);
        auto dec1 = decrypt_scalar(sk1, ct1s, params);
        bool ok = true;
        auto dec2 = decrypt_vector(sk2, ct2s, params, ok);
        bool s1 = (dec1 == mu1);
        bool s2 = ok && (dec2[0] == m2[0]) && (dec2[1] == m2[1]);
        cout << "[sanity] scalar dec: " << (s1 ? "OK" : "FAIL")
             << ",  vector dec: " << (s2 ? "OK" : "FAIL")
             << ",  sparsification check: " << (ok ? "passed" : "FAILED") << "\n";
    }

    // ---------- Multi-key Mul at full level ----------
    auto t_mul0 = std::chrono::steady_clock::now();
    auto tct = mkMul(ct1, ct2, params);
    auto t_mul1 = std::chrono::steady_clock::now();
    cout << "[mkMul]  tensor ct components: 1 (d0) + " << tct.d1.size()
         << " (d1) + 1 (d2) + " << tct.d3.size() << " (d3), L=" << tct.d0.num_primes
         << ",  " << std::chrono::duration_cast<std::chrono::milliseconds>(t_mul1 - t_mul0).count() << " ms\n";

    // ---------- Manual mod switch: L=3 -> L=2 -> L=1  (log q' ~ 44 bits) ----------
    auto t_ms0 = std::chrono::steady_clock::now();
    mod_switch(tct, params);
    cout << "[modswitch 1] now L=" << tct.d0.num_primes
         << ", remaining log q = "
         << (bit_length(params.q_chain[0].value()) + bit_length(params.q_chain[1].value())) << " bits\n";
    mod_switch(tct, params);
    cout << "[modswitch 2] now L=" << tct.d0.num_primes
         << ", log q' = " << bit_length(params.q_chain[0].value()) << " bits  (target ~44)\n";
    auto t_ms1 = std::chrono::steady_clock::now();
    cout << "[modswitch] " << std::chrono::duration_cast<std::chrono::milliseconds>(t_ms1 - t_ms0).count() << " ms\n";

    // ---------- Multi-key Decryption ----------
    auto t_dec0 = std::chrono::steady_clock::now();
    bool mk_check = false;
    auto mk_out = mkDec(sk1, sk2, tct, params, mk_check);
    auto t_dec1 = std::chrono::steady_clock::now();
    cout << "[mkDec]  " << std::chrono::duration_cast<std::chrono::milliseconds>(t_dec1 - t_dec0).count() << " ms,"
         << "  sparsification check: " << (mk_check ? "passed" : "FAILED") << "\n";

    // ---------- Verification ----------
    bool all_ok = mk_check;
    for (size_t c = 0; c < params.ell; c++) {
        auto expected = plain_mul(mu1, m2[c], params);
        bool match = (expected == mk_out[c]);
        cout << "[verify] component " << c << " (= μ1 * m2[" << c << "] in R_p):  "
             << (match ? "MATCH" : "MISMATCH") << "\n";
        if (!match) {
            size_t mismatches = 0;
            for (size_t i = 0; i < params.N; i++) if (expected[i] != mk_out[c][i]) mismatches++;
            cout << "         coefficient mismatches: " << mismatches << " / " << params.N << "\n";
            cout << "         expected[0..3]: ";
            for (int i = 0; i < 4; i++) cout << expected[i] << " ";
            cout << "\n         got     [0..3]: ";
            for (int i = 0; i < 4; i++) cout << mk_out[c][i] << " ";
            cout << "\n";
            all_ok = false;
        }
    }

    cout << "\n========== " << (all_ok ? "SUCCESS" : "FAILURE") << " ==========\n";
    return all_ok ? 0 : 1;
}
