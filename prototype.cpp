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

} // namespace mk

// ============================== main ==============================
int main() {
    using namespace mk;
    using std::cout;
    using std::endl;

    cout << "=== Multi-key RLWE Vector Encryption — BGV prototype (one mul layer) ===\n";

    Params params;
    const uint64_t plain_modulus = (uint64_t(7) << 20) + 1; // 7340033 = 7*2^20 + 1
    // Three primes ≈ {44, 43, 43} bits, total ≈ 130 bits; manual mod-switch twice -> 1 prime ≈ 44 bits.
    vector<int> q_bits = {44, 43, 43};
    params_init(params, plain_modulus, q_bits, /*N*/ 8192, /*ell*/ 2, /*tau*/ 1, /*hw*/ 32);

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
