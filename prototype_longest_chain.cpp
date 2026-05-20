// LONGEST-CHAIN sweep test (non-ZK).
//
// Parameterized version: N_OPS = total number of BGV operations before mkMul.
// Structure of the N_OPS BGV ops:
//   - op 1:           F_p × ct1
//   - op 2:           ct1 × ct1 + relin
//   - ops 3..N_OPS:   F_p × ct1   (plain multiplications)
// Mod-switch every 2 ops. No smudge (no ZK).
// Final tensor mod-switched L=3 → L=1 just like the other tests.
//
// Set N_OPS via env var. Defaults to 4 (matches prototype_full_chain_nozk).
// Use SWEEP=1 to run a sweep over N_OPS ∈ {4, 6, 8, 10, 14, 20, 30, 50, 100} and report
// initial budget + decrypt success for each.

#include <algorithm>
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

#include "seal/seal.h"
#include "seal/util/numth.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/ntt.h"
#include "seal/util/uintarithsmallmod.h"

using namespace seal;
using namespace seal::util;
using namespace std;
using namespace std::chrono;

// ============================================================================
// Prototype framework (minimal copy from prototype.cpp)
// ============================================================================
namespace mk {

class Prng {
public:
    Prng() : rng_(random_device{}()) {}
    explicit Prng(uint64_t seed) : rng_(seed) {}
    uint64_t u64() { return rng_(); }
    uint64_t uniform_below(uint64_t bound) { return bound ? (rng_() % bound) : 0; }
    int cbd(int eta) {
        int s = 0;
        for (int i = 0; i < eta; i++) s += int(rng_() & 1ULL);
        for (int i = 0; i < eta; i++) s -= int(rng_() & 1ULL);
        return s;
    }
private:
    mt19937_64 rng_;
};

struct Params {
    size_t N = 0;
    int log_N = 0;
    Modulus p;
    vector<Modulus> q_chain;
    vector<unique_ptr<NTTTables>> q_ntt;
    unique_ptr<NTTTables> p_ntt;
    size_t ell = 0, tau = 0, ell_prime = 0;
    size_t hw = 0;
    vector<uint64_t> p_mod_q;
    int sigma_eta = 21;
};

static int bit_length(uint64_t v) { int b = 0; while (v) { b++; v >>= 1; } return b; }

static void params_init_from_moduli(Params& params, uint64_t plain_modulus_val,
                                     const vector<Modulus>& moduli, size_t N,
                                     size_t ell, size_t tau, size_t hw) {
    params.N = N;
    params.log_N = 0;
    while ((size_t(1) << params.log_N) < N) params.log_N++;
    params.p = Modulus(plain_modulus_val);
    params.ell = ell;
    params.tau = tau;
    params.ell_prime = ell + tau;
    params.hw = hw;
    params.q_chain = moduli;
    params.q_ntt.clear();
    for (auto& q : params.q_chain) {
        params.q_ntt.emplace_back(make_unique<NTTTables>(params.log_N, q));
    }
    params.p_ntt = make_unique<NTTTables>(params.log_N, params.p);
    params.p_mod_q.resize(params.q_chain.size());
    for (size_t i = 0; i < params.q_chain.size(); i++) {
        params.p_mod_q[i] = barrett_reduce_64(plain_modulus_val, params.q_chain[i]);
    }
}

struct RnsPoly {
    vector<uint64_t> data;
    size_t num_primes = 0;
    size_t N = 0;
};

static inline uint64_t* row(RnsPoly& p, size_t k) { return p.data.data() + k * p.N; }
static inline const uint64_t* row(const RnsPoly& p, size_t k) { return p.data.data() + k * p.N; }

static RnsPoly make_zero(const Params& params, size_t L) {
    RnsPoly r; r.num_primes = L; r.N = params.N; r.data.assign(L * params.N, 0ULL); return r;
}

static RnsPoly sample_uniform(const Params& params, size_t L, Prng& prng) {
    RnsPoly r = make_zero(params, L);
    for (size_t k = 0; k < L; k++) {
        uint64_t q = params.q_chain[k].value();
        uint64_t* dst = row(r, k);
        for (size_t i = 0; i < params.N; i++) dst[i] = prng.u64() % q;
    }
    return r;
}

static vector<int> sample_ternary_hw(size_t N, size_t hw, Prng& prng) {
    vector<int> coeffs(N, 0);
    vector<size_t> positions(N);
    for (size_t i = 0; i < N; i++) positions[i] = i;
    for (size_t i = 0; i < hw; i++) {
        size_t j = i + prng.uniform_below(N - i);
        swap(positions[i], positions[j]);
    }
    for (size_t i = 0; i < hw; i++) coeffs[positions[i]] = (prng.u64() & 1ULL) ? 1 : -1;
    return coeffs;
}

static RnsPoly small_signed_to_rns_ntt(const vector<int>& coeffs, const Params& params, size_t L) {
    RnsPoly r = make_zero(params, L);
    for (size_t k = 0; k < L; k++) {
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

static RnsPoly plain_to_rns_ntt(const vector<uint64_t>& coeffs, const Params& params, size_t L) {
    RnsPoly r = make_zero(params, L);
    for (size_t k = 0; k < L; k++) {
        uint64_t* dst = row(r, k);
        memcpy(dst, coeffs.data(), params.N * sizeof(uint64_t));
        ntt_negacyclic_harvey(dst, *params.q_ntt[k]);
    }
    return r;
}

static vector<int> sample_noise_cbd(size_t N, int eta, Prng& prng) {
    vector<int> coeffs(N);
    for (size_t i = 0; i < N; i++) coeffs[i] = prng.cbd(eta);
    return coeffs;
}

static RnsPoly add(const RnsPoly& a, const RnsPoly& b, const Params& params) {
    RnsPoly r = make_zero(params, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++)
        add_poly_coeffmod(row(a, k), row(b, k), a.N, params.q_chain[k], row(r, k));
    return r;
}
static RnsPoly sub(const RnsPoly& a, const RnsPoly& b, const Params& params) {
    RnsPoly r = make_zero(params, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++)
        sub_poly_coeffmod(row(a, k), row(b, k), a.N, params.q_chain[k], row(r, k));
    return r;
}
static RnsPoly mul(const RnsPoly& a, const RnsPoly& b, const Params& params) {
    RnsPoly r = make_zero(params, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++)
        dyadic_product_coeffmod(row(a, k), row(b, k), a.N, params.q_chain[k], row(r, k));
    return r;
}
static RnsPoly mul_by_p(const RnsPoly& a, const Params& params) {
    RnsPoly r = make_zero(params, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++) {
        multiply_poly_scalar_coeffmod(row(a, k), a.N, params.p_mod_q[k], params.q_chain[k], row(r, k));
    }
    return r;
}

static RnsPoly truncate_to_level(const RnsPoly& s, const Params& params, size_t target_primes) {
    RnsPoly r;
    r.N = params.N;
    r.num_primes = target_primes;
    r.data.assign(s.data.begin(), s.data.begin() + target_primes * params.N);
    return r;
}

static void mod_switch_drop_last_poly(RnsPoly& c, const Params& params) {
    if (c.num_primes <= 1) throw runtime_error("cannot drop last prime: only one left");
    const size_t L = c.num_primes;
    const Modulus& q_last = params.q_chain[L - 1];
    const uint64_t q_last_val = q_last.value();
    vector<uint64_t> c_last(params.N);
    memcpy(c_last.data(), row(c, L - 1), params.N * sizeof(uint64_t));
    inverse_ntt_negacyclic_harvey(c_last.data(), *params.q_ntt[L - 1]);
    vector<uint64_t> ncl_mod_p(params.N);
    modulo_poly_coeffs(c_last.data(), params.N, params.p, ncl_mod_p.data());
    negate_poly_coeffmod(ncl_mod_p.data(), params.N, params.p, ncl_mod_p.data());
    uint64_t inv_q_last_mod_p = 1;
    if (!try_invert_uint_mod(q_last_val, params.p, inv_q_last_mod_p))
        throw runtime_error("q_last not invertible mod p");
    if (inv_q_last_mod_p != 1) {
        multiply_poly_scalar_coeffmod(ncl_mod_p.data(), params.N, inv_q_last_mod_p, params.p, ncl_mod_p.data());
    }
    vector<uint64_t> delta(params.N);
    for (size_t k = 0; k < L - 1; k++) {
        const Modulus& q_k = params.q_chain[k];
        modulo_poly_coeffs(ncl_mod_p.data(), params.N, q_k, delta.data());
        uint64_t qlast_mod_qk = barrett_reduce_64(q_last_val, q_k);
        multiply_poly_scalar_coeffmod(delta.data(), params.N, qlast_mod_qk, q_k, delta.data());
        for (size_t i = 0; i < params.N; i++) {
            uint64_t cl_mod_qk = barrett_reduce_64(c_last[i], q_k);
            delta[i] = add_uint_mod(delta[i], cl_mod_qk, q_k);
        }
        ntt_negacyclic_harvey(delta.data(), *params.q_ntt[k]);
        sub_poly_coeffmod(row(c, k), delta.data(), params.N, q_k, row(c, k));
        uint64_t inv_qlast_mod_qk = 0;
        if (!try_invert_uint_mod(q_last_val, q_k, inv_qlast_mod_qk))
            throw runtime_error("q_last not invertible mod q_k");
        multiply_poly_scalar_coeffmod(row(c, k), params.N, inv_qlast_mod_qk, q_k, row(c, k));
    }
    c.data.resize((L - 1) * params.N);
    c.num_primes = L - 1;
}

static vector<uint64_t> reduce_to_plain(const RnsPoly& c, const Params& params, uint64_t correction) {
    const size_t L = c.num_primes;
    vector<uint64_t> out(params.N);
    if (L == 1) {
        vector<uint64_t> coeffs(params.N);
        memcpy(coeffs.data(), c.data.data(), params.N * sizeof(uint64_t));
        inverse_ntt_negacyclic_harvey(coeffs.data(), *params.q_ntt[0]);
        const uint64_t q = params.q_chain[0].value();
        const uint64_t q_half = q >> 1;
        const uint64_t p_val = params.p.value();
        const uint64_t q_mod_p = barrett_reduce_64(q, params.p);
        for (size_t i = 0; i < params.N; i++) {
            uint64_t v = coeffs[i];
            uint64_t v_mod_p = barrett_reduce_64(v, params.p);
            if (v > q_half) {
                if (v_mod_p >= q_mod_p) out[i] = v_mod_p - q_mod_p;
                else out[i] = p_val - (q_mod_p - v_mod_p);
            } else {
                out[i] = v_mod_p;
            }
        }
    } else if (L == 2) {
        vector<uint64_t> c0(params.N), c1(params.N);
        memcpy(c0.data(), row(c, 0), params.N * sizeof(uint64_t));
        memcpy(c1.data(), row(c, 1), params.N * sizeof(uint64_t));
        inverse_ntt_negacyclic_harvey(c0.data(), *params.q_ntt[0]);
        inverse_ntt_negacyclic_harvey(c1.data(), *params.q_ntt[1]);
        const uint64_t q0 = params.q_chain[0].value();
        const uint64_t q1 = params.q_chain[1].value();
        const __uint128_t Q = __uint128_t(q0) * __uint128_t(q1);
        const __uint128_t Q_half = Q >> 1;
        uint64_t inv_q0_mod_q1;
        if (!try_invert_uint_mod(q0, params.q_chain[1], inv_q0_mod_q1))
            throw runtime_error("q0 not invertible mod q1");
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
        throw runtime_error("reduce_to_plain: only L in {1,2} supported");
    }
    if (correction != 1) {
        uint64_t inv_corr;
        if (!try_invert_uint_mod(correction, params.p, inv_corr))
            throw runtime_error("correction not invertible mod p");
        multiply_poly_scalar_coeffmod(out.data(), params.N, inv_corr, params.p, out.data());
    }
    return out;
}

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

// ---- key & ct types ----
struct SecretKeyScalar { vector<int> s_coeffs; RnsPoly s_ntt; };
struct SecretKeyVector {
    vector<vector<int>> s_coeffs;
    vector<RnsPoly> s_ntt;
    vector<vector<vector<uint64_t>>> T_coeffs;  // R_p^(tau × ell)
};
struct CipherScalar { RnsPoly a, b; uint64_t correction = 1; };
struct CipherVector { RnsPoly a; vector<RnsPoly> b; uint64_t correction = 1; };
struct TensorCt {
    RnsPoly d0;
    vector<RnsPoly> d1;
    RnsPoly d2;
    vector<RnsPoly> d3;
    uint64_t correction = 1;
};

static SecretKeyVector keygen_vector(const Params& params, Prng& prng) {
    SecretKeyVector sk;
    const size_t L = params.q_chain.size();
    sk.s_coeffs.resize(params.ell_prime);
    sk.s_ntt.reserve(params.ell_prime);
    for (size_t i = 0; i < params.ell_prime; i++) {
        sk.s_coeffs[i] = sample_ternary_hw(params.N, params.hw, prng);
        sk.s_ntt.emplace_back(small_signed_to_rns_ntt(sk.s_coeffs[i], params, L));
    }
    sk.T_coeffs.assign(params.tau, vector<vector<uint64_t>>(params.ell, vector<uint64_t>(params.N, 0)));
    for (size_t r = 0; r < params.tau; r++) {
        for (size_t c = 0; c < params.ell; c++) {
            for (size_t i = 0; i < params.N; i++) sk.T_coeffs[r][c][i] = prng.uniform_below(params.p.value());
        }
    }
    return sk;
}

static CipherVector encrypt_vector(const SecretKeyVector& sk,
                                    const vector<vector<uint64_t>>& msgs_m,
                                    const Params& params, Prng& prng) {
    const size_t L = params.q_chain.size();
    vector<vector<uint64_t>> mprime(params.tau, vector<uint64_t>(params.N, 0));
    for (size_t r = 0; r < params.tau; r++) {
        for (size_t c = 0; c < params.ell; c++) {
            auto prod = plain_mul(sk.T_coeffs[r][c], msgs_m[c], params);
            mprime[r] = plain_add(mprime[r], prod, params);
        }
    }
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

static TensorCt mkMul(const CipherScalar& ct1, const CipherVector& ct2, const Params& params) {
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

static void mod_switch_tensor(TensorCt& ct, const Params& params) {
    if (ct.d0.num_primes <= 1) throw runtime_error("cannot mod-switch tensor");
    const Modulus& q_last = params.q_chain[ct.d0.num_primes - 1];
    uint64_t corr;
    if (!try_invert_uint_mod(q_last.value(), params.p, corr))
        throw runtime_error("q_last not invertible mod p");
    mod_switch_drop_last_poly(ct.d0, params);
    mod_switch_drop_last_poly(ct.d2, params);
    for (auto& d : ct.d1) mod_switch_drop_last_poly(d, params);
    for (auto& d : ct.d3) mod_switch_drop_last_poly(d, params);
    ct.correction = multiply_uint_mod(ct.correction, corr, params.p);
}

static vector<vector<uint64_t>> mkDec(const SecretKeyScalar& sk1, const SecretKeyVector& sk2,
                                      const TensorCt& tct, const Params& params,
                                      bool& sparsi_ok) {
    const size_t L = tct.d0.num_primes;
    RnsPoly s1_lvl = truncate_to_level(sk1.s_ntt, params, L);
    vector<vector<uint64_t>> mu_full(params.ell_prime);
    for (size_t i = 0; i < params.ell_prime; i++) {
        RnsPoly s2i_lvl = truncate_to_level(sk2.s_ntt[i], params, L);
        RnsPoly t1 = mul(mul(tct.d0, s1_lvl, params), s2i_lvl, params);
        RnsPoly t2 = mul(tct.d1[i], s1_lvl, params);
        RnsPoly t3 = mul(tct.d2, s2i_lvl, params);
        RnsPoly r = add(sub(sub(t1, t2, params), t3, params), tct.d3[i], params);
        mu_full[i] = reduce_to_plain(r, params, tct.correction);
    }
    sparsi_ok = true;
    for (size_t r = 0; r < params.tau; r++) {
        vector<uint64_t> expected(params.N, 0);
        for (size_t c = 0; c < params.ell; c++) {
            auto prod = plain_mul(sk2.T_coeffs[r][c], mu_full[c], params);
            expected = plain_add(expected, prod, params);
        }
        if (expected != mu_full[params.ell + r]) sparsi_ok = false;
    }
    vector<vector<uint64_t>> out(params.ell);
    for (size_t i = 0; i < params.ell; i++) out[i] = mu_full[i];
    return out;
}

static void add_smudge_poly(RnsPoly& target, int smudge_log2, const Params& params, Prng& prng) {
    if (smudge_log2 < 1 || smudge_log2 > 62) throw runtime_error("smudge_log2 out of range");
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

// Add p · Gaussian(σ) noise to a tensor component, simulating noise from sum of K fresh outputs.
static void add_gaussian_p_noise(RnsPoly& target, double sigma, const Params& params, Prng& prng) {
    const size_t L = target.num_primes;
    std::mt19937_64 rng(prng.u64());
    std::normal_distribution<double> nd(0.0, sigma);
    vector<int64_t> n_signed(params.N);
    for (size_t i = 0; i < params.N; i++) n_signed[i] = (int64_t)std::llround(nd(rng));
    RnsPoly nz = make_zero(params, L);
    for (size_t k = 0; k < L; k++) {
        uint64_t q = params.q_chain[k].value();
        uint64_t* dst = row(nz, k);
        for (size_t i = 0; i < params.N; i++) {
            int64_t v = n_signed[i];
            if (v < 0) {
                int64_t r = (-v) % (int64_t)q;
                dst[i] = (r == 0) ? 0 : (q - (uint64_t)r);
            } else {
                dst[i] = (uint64_t)(v % (int64_t)q);
            }
        }
        ntt_negacyclic_harvey(dst, *params.q_ntt[k]);
    }
    RnsPoly p_n = mul_by_p(nz, params);
    target = add(target, p_n, params);
}

} // namespace mk

// ============================================================================
// SEAL → prototype bridge
// ============================================================================
// Extract SEAL's secret key into our prototype's RnsPoly form (at L=3, the first 3 prime slabs).
static mk::SecretKeyScalar extract_sk_to_proto(const SecretKey& sk_seal,
                                                const SEALContext& ctx_seal,
                                                const mk::Params& params_proto) {
    const size_t N = params_proto.N;
    const size_t L = params_proto.q_chain.size();
    mk::SecretKeyScalar sk;
    sk.s_ntt.N = N;
    sk.s_ntt.num_primes = L;
    sk.s_ntt.data.resize(L * N);
    // sk.data() is a Plaintext storing s in NTT form against the KEY context.
    // Key context primes = data primes + special prime. Order: data primes first, special last.
    // So bottom L primes of key context = primes of L=3 context. Copy L * N uint64s from the front.
    memcpy(sk.s_ntt.data.data(), sk_seal.data().data(), L * N * sizeof(uint64_t));
    (void)ctx_seal;
    return sk;
}

// Convert SEAL Ciphertext (after mod-switch chain) into our prototype's CipherScalar.
// Computes the correction = product of (q_dropped^{-1} mod p) over all primes dropped from
// SEAL's first_parms_id chain down to the ct's current level.
static mk::CipherScalar seal_ct_to_proto(const Ciphertext& seal_ct,
                                          const SEALContext& ctx_seal,
                                          const mk::Params& params_proto) {
    const size_t N = params_proto.N;
    const size_t L = params_proto.q_chain.size();
    if (seal_ct.coeff_modulus_size() != L) {
        throw runtime_error("seal_ct_to_proto: ct level doesn't match prototype params");
    }
    if (!seal_ct.is_ntt_form()) {
        throw runtime_error("seal_ct_to_proto: SEAL BGV ct expected in NTT form");
    }

    mk::CipherScalar out;
    out.a = mk::make_zero(params_proto, L);
    out.b = mk::make_zero(params_proto, L);

    // Convention bridge: SEAL c0 = b_seal = -a·s + p·e + μ, c1 = a_seal.
    // Prototype: b_proto = a_proto·s + p·e + μ. Set a_proto = −a_seal (negate c1), b_proto = c0.
    for (size_t k = 0; k < L; k++) {
        uint64_t q = params_proto.q_chain[k].value();
        uint64_t* dst_a = mk::row(out.a, k);
        uint64_t* dst_b = mk::row(out.b, k);
        const uint64_t* src_c0 = seal_ct.data(0) + k * N;
        const uint64_t* src_c1 = seal_ct.data(1) + k * N;
        for (size_t i = 0; i < N; i++) {
            dst_a[i] = (src_c1[i] == 0) ? 0 : (q - src_c1[i]);  // negate mod q
            dst_b[i] = src_c0[i];
        }
    }

    // SEAL tracks the BGV correction factor in Ciphertext::correction_factor(), keeping it
    // synchronized through mod-switches/mul. Use it directly.
    out.correction = seal_ct.correction_factor();
    (void)ctx_seal;
    return out;
}

// ============================================================================
// Main pipeline
// ============================================================================
int main() {
    const size_t N = 8192;
    const uint64_t p_val = (uint64_t(7) << 20) + 1;          // 7340033

    // SEAL chain (bottom→top, SEAL ordering with special prime last):
    //   q0=33  (final prime after step 10)
    //   q1=49, q2=48 (dropped during step 10)
    //   q3..q7=43 each (dropped during steps 1-5 intermediate mod-switches)
    //   q_special=60
    // Top-level data q ≈ 345 bits → initial noise budget ~318 bits.
    EncryptionParameters parms(scheme_type::bgv);
    parms.set_poly_modulus_degree(N);
    // Default final prime 34b (min needed for the chain). Override via Q0_BITS.
    int q0_bits = 34;
    if (const char* e = getenv("Q0_BITS")) q0_bits = atoi(e);
    int q1_bits = 49, q2_bits = 47;
    if (q0_bits != 34) {
        q2_bits = 130 - q0_bits - q1_bits;
        if (q2_bits < 33) q2_bits = 33;
    }
    // N_OPS env var = total BGV operations (op 1 plain + op 2 ct×ct + (N_OPS-2) plain muls).
    // Must be even. Defaults to 4. With mod-switch every 2 ops, L_top = N_OPS/2 + 3.
    int n_ops = 4;
    if (const char* e = getenv("N_OPS")) n_ops = atoi(e);
    if (n_ops < 2) n_ops = 2;
    if (n_ops % 2 != 0) n_ops += 1;     // round up to even

    const int n_upper_primes = n_ops / 2;   // number of 49-bit upper primes
    const int L_top = n_upper_primes + 3;   // data primes total

    // Build coeff_modulus: {q0, q1, q2, 49, 49, ..., 49, 60_special}
    vector<int> chain_bits = {(int)q0_bits, q1_bits, q2_bits};
    for (int i = 0; i < n_upper_primes; i++) chain_bits.push_back(49);
    chain_bits.push_back(60);   // special prime
    parms.set_coeff_modulus(CoeffModulus::Create(N, chain_bits));
    parms.set_plain_modulus(p_val);
    SEALContext ctx(parms, /*expand_mod_chain*/ true, sec_level_type::none);

    KeyGenerator kg(ctx);
    SecretKey sk = kg.secret_key();
    PublicKey pk; kg.create_public_key(pk);
    RelinKeys rk; kg.create_relin_keys(rk);
    Encryptor enc(ctx, pk);
    enc.set_secret_key(sk);
    Decryptor dec(ctx, sk);
    Evaluator eval(ctx);
    BatchEncoder batcher(ctx);

    // Build prototype params for the L=3 context (3 primes, matching ct2).
    // We'll fill this in after we walk the chain down to L=3.
    mk::Params proto_params;
    {
        auto cd = ctx.first_context_data();
        // Walk down to L=3: drop 5 primes from the data chain.
        size_t target_primes = 3;
        while (cd->parms().coeff_modulus().size() > target_primes) cd = cd->next_context_data();
        const auto& mods = cd->parms().coeff_modulus();
        mk::params_init_from_moduli(proto_params, p_val, mods, N,
                                      /*ell*/ 4, /*tau*/ 1, /*hw*/ 32);
    }

    // Report params.
    cout << "=== Full chain: BGV → multi-key → smudge → mod-switch ===\n";
    cout << "  N = " << N << ",  p = " << p_val << "  (~" << mk::bit_length(p_val) << " b)\n";
    cout << "  SEAL chain (data primes + special):  ";
    int top_bits = 0;
    auto cd_top = ctx.first_context_data();
    for (auto& m : cd_top->parms().coeff_modulus()) {
        cout << mk::bit_length(m.value()) << "b ";
        top_bits += mk::bit_length(m.value());
    }
    int sp_bits = 0;
    {
        auto& key_mods = ctx.key_context_data()->parms().coeff_modulus();
        sp_bits = mk::bit_length(key_mods.back().value());
    }
    cout << "+ " << sp_bits << "b(special)";
    cout << "\n  Top-level data q  ≈ " << top_bits << " bits\n";
    cout << "  L=3 q (after step 6) ≈ ";
    {
        int l3_bits = 0;
        for (auto& m : proto_params.q_chain) l3_bits += mk::bit_length(m.value());
        cout << l3_bits << " bits  ({";
        for (size_t k = 0; k < proto_params.q_chain.size(); k++)
            cout << (k ? ", " : "") << mk::bit_length(proto_params.q_chain[k].value()) << "b";
        cout << "})\n";
    }
    cout << "\n";

    // ---- Encrypt ct1 with μ = constant 1 polynomial (coefficient form: [1, 0, 0, ..., 0]) ----
    // We deliberately DON'T use BatchEncoder, since it applies a slot transform: encoding
    // [1, 1, ..., 1] would give a Plaintext whose polynomial is constant 1, but the prototype's
    // reduce_to_plain works in coefficient form, so it'd give [1, 0, ..., 0] back — confusing.
    // Direct raw Plaintext keeps polynomial = coefficient throughout.
    Plaintext mu1_pt("1");   // polynomial: 1 (= 1 + 0·X + 0·X^2 + ...)
    Ciphertext ct1;
    enc.encrypt_symmetric(mu1_pt, ct1);

    int budget_init = dec.invariant_noise_budget(ct1);
    cout << "[init] Initial noise budget for ct_1: " << budget_init << " bits\n\n";

    // For plain-mul steps, use raw "1" polynomial too — so ct1's polynomial stays = 1
    // throughout the chain (1 × 1 = 1, 1^2 = 1).
    Plaintext m_pt("1");
    eval.transform_to_ntt_inplace(m_pt, ct1.parms_id());

    auto t_start = steady_clock::now();
    auto take_us = [&]() {
        auto now = steady_clock::now();
        auto us = duration_cast<microseconds>(now - t_start).count();
        t_start = now;
        return (long)us;
    };
    auto fmt = [](long us) {
        ostringstream o; o << setw(8) << us;
        return o.str();
    };

    cout << "Stage                                     us       budget(b)\n";
    cout << "-----                                  ------     ----------\n";

    // Parameterized chain with N_OPS BGV operations:
    //   op 1: F_p × ct1
    //   op 2: ct1 × ct1 + relin
    //   ops 3..N_OPS: F_p × ct1
    // Mod-switch happens after EVERY 2nd op. No smudge applied to the tensor (non-ZK).
    cout << "Running with N_OPS = " << n_ops
         << ",  L_top = " << L_top << " data primes,  upper primes: "
         << n_upper_primes << " × 49b\n\n";

    // op 1: F_p × ct1
    t_start = steady_clock::now();
    eval.multiply_plain_inplace(ct1, m_pt);
    long t1 = take_us();
    int b1 = dec.invariant_noise_budget(ct1);
    cout << "  op 1  (F_p × ct_1)                    " << fmt(t1) << "       " << b1 << "\n";

    // op 2: ct1 × ct1 + relin
    Ciphertext ct1_copy = ct1;
    t_start = steady_clock::now();
    eval.multiply_inplace(ct1, ct1_copy);
    long t2 = take_us();
    int b2 = dec.invariant_noise_budget(ct1);
    cout << "  op 2a (ct_1 × ct_1)                   " << fmt(t2) << "       " << b2 << "\n";

    t_start = steady_clock::now();
    eval.relinearize_inplace(ct1, rk);
    long t2r = take_us();
    int b2r = dec.invariant_noise_budget(ct1);
    cout << "  op 2b (relinearize)                   " << fmt(t2r) << "       " << b2r << "\n";

    // mod-switch after the first pair (ops 1, 2)
    t_start = steady_clock::now();
    eval.mod_switch_to_next_inplace(ct1);
    long t2m = take_us();
    int b2m = dec.invariant_noise_budget(ct1);
    cout << "    + mod-switch (after pair 1)         " << fmt(t2m) << "       " << b2m << "\n";

    // ops 3..N_OPS: plain muls, mod-switch every 2 ops
    for (int op = 3; op <= n_ops; op++) {
        Plaintext mp("1");
        eval.transform_to_ntt_inplace(mp, ct1.parms_id());
        t_start = steady_clock::now();
        eval.multiply_plain_inplace(ct1, mp);
        long t = take_us();
        int b = dec.invariant_noise_budget(ct1);
        cout << "  op " << setw(2) << op << " (F_p × ct_1)                   " << fmt(t) << "       " << b << "\n";

        // Mod-switch after every 2 ops (i.e. after even-numbered ops).
        if (op % 2 == 0) {
            t_start = steady_clock::now();
            eval.mod_switch_to_next_inplace(ct1);
            long tm = take_us();
            int bm = dec.invariant_noise_budget(ct1);
            cout << "    + mod-switch (after pair " << (op / 2) << ")         "
                 << fmt(tm) << "       " << bm << "\n";
        }
    }

    // Drain any remaining levels to reach L=3 if not already there (handles odd N_OPS edge cases).
    {
        size_t cur = ct1.coeff_modulus_size();
        cout << "  Step 6 (drain to L=3)";
        if (cur > 3) {
            t_start = steady_clock::now();
            while (ct1.coeff_modulus_size() > 3) eval.mod_switch_to_next_inplace(ct1);
            long t = take_us();
            int b = dec.invariant_noise_budget(ct1);
            cout << "                  " << fmt(t) << "       " << b << "\n";
        } else {
            cout << "                  (no-op)   " << dec.invariant_noise_budget(ct1) << "\n";
        }
    }

    // --- Verify SEAL ct1 still decrypts correctly ---
    // Underlying polynomial is constant 1: pt_check.coeff_count()==1 and pt_check[0]==1.
    {
        Plaintext pt_check; dec.decrypt(ct1, pt_check);
        // Constant 1 polynomial: pt_check should equal "1"
        bool ok = (pt_check.to_string() == string("1"));
        cout << "  [post-step-6 SEAL decrypt check]   " << (ok ? "OK" : "FAIL")
             << "  (expected polynomial = constant 1; got \"" << pt_check.to_string() << "\")\n\n";
    }

    // ---- Bridge: SEAL ct1 → prototype CipherScalar ----
    mk::CipherScalar ct1_proto = seal_ct_to_proto(ct1, ctx, proto_params);
    mk::SecretKeyScalar sk1_proto = extract_sk_to_proto(sk, ctx, proto_params);

    // Debug: copy ct1_proto, mod-switch L=3 → L=2 in the prototype framework, decrypt at L=2.
    // Verifies bridge + sk extraction + mod-switch all work for the scalar case.
    {
        mk::CipherScalar ct1_test = ct1_proto;
        // Inline mod_switch for CipherScalar: drop the last prime from a and b, update correction.
        const Modulus& q_last = proto_params.q_chain[ct1_test.a.num_primes - 1];
        uint64_t corr;
        try_invert_uint_mod(q_last.value(), proto_params.p, corr);
        mk::mod_switch_drop_last_poly(ct1_test.a, proto_params);
        mk::mod_switch_drop_last_poly(ct1_test.b, proto_params);
        ct1_test.correction = multiply_uint_mod(ct1_test.correction, corr, proto_params.p);

        mk::RnsPoly s_lvl = mk::truncate_to_level(sk1_proto.s_ntt, proto_params, ct1_test.a.num_primes);
        mk::RnsPoly as = mk::mul(ct1_test.a, s_lvl, proto_params);
        mk::RnsPoly bma = mk::sub(ct1_test.b, as, proto_params);
        auto dec_msg = mk::reduce_to_plain(bma, proto_params, ct1_test.correction);
        bool ok = true;
        size_t first_bad = N;
        for (size_t i = 0; i < N; i++) {
            uint64_t expected = (i == 0) ? 1 : 0;   // polynomial = constant 1 → coeff form [1,0,0,...]
            if (dec_msg[i] != expected) { ok = false; first_bad = i; break; }
        }
        cout << "  [debug: prototype mod-switch + decrypt ct1_proto at L=2] "
             << (ok ? "OK" : "FAIL");
        if (!ok) cout << " (first bad: idx " << first_bad << ", got " << dec_msg[first_bad] << ")";
        cout << "\n";
    }

    // ---- Setup ct2 in prototype framework ----
    mk::Prng prng_proto(20260519ULL);
    auto sk2 = mk::keygen_vector(proto_params, prng_proto);
    vector<vector<uint64_t>> m2_data(proto_params.ell, vector<uint64_t>(N));
    for (size_t i = 0; i < proto_params.ell; i++)
        for (size_t j = 0; j < N; j++) m2_data[i][j] = prng_proto.uniform_below(p_val);
    auto ct2_proto = mk::encrypt_vector(sk2, m2_data, proto_params, prng_proto);

    // (Skip ct2 sanity decrypt at L=3 for the same reason.)
    cout << "\n";

    // ============================================================
    // Steps 7-10 in prototype framework
    // ============================================================
    cout << "Stage                                     us       log2||e||\n";
    cout << "-----                                  ------     ----------\n";

    // Step 7: mkMul
    t_start = steady_clock::now();
    auto tct = mk::mkMul(ct1_proto, ct2_proto, proto_params);
    long t7 = take_us();
    cout << "  Step 7 (mkMul → 12-poly tensor)       " << fmt(t7) << "\n";

    // Step 8: simulate sum of K=2^20 independent fresh ciphertexts (same μ).
    //   - Plaintext should scale by K: achieved by ADDING (K-1)·μ·correction to d3[i] (lifted to R_q
    //     in NTT form). This contributes to the plaintext side of mkDec, not the noise side.
    //   - Noise grows as random walk √K in std-dev: add Gaussian noise with σ = √(K-1)·σ_per to d3.
    //     d3 enters mkDec without any s-multiplication (just `+ d3`), so noise lands cleanly.
    //   - σ_per (per-coefficient std of fresh-mkMul noise after a deep mod-switch): empirically
    //     observed to be small after mod-switch absorbs most of the noise. We use ~2^32 as a
    //     safe upper-bound estimate; the result is checked by mkDec at end.
    const uint64_t K = (uint64_t(1) << 20);
    const uint64_t K_minus_1_mod_p = (K - 1) % p_val;
    const double sigma_per = std::ldexp(1.0, 32);             // ≈ 2^32 (rough upper bound)
    const double sigma_K   = std::sqrt((double)(K - 1)) * sigma_per;  // ≈ 2^42

    t_start = steady_clock::now();
    // [DEBUG] To isolate the bridge, temporarily SKIP step 8's plaintext scaling and noise add.
    // We'll re-enable once basic decryption is verified.
    const bool SKIP_STEP_8 = (getenv("SKIP_STEP_8") != nullptr);
    const bool SKIP_STEP_9 = (getenv("SKIP_STEP_9") != nullptr);
    for (size_t i = 0; i < proto_params.ell_prime; i++) {
        if (SKIP_STEP_8) break;
        // Compute the underlying μ_i for slot i (= m2_data[i] for i<ell, else T·m for sparsi rows).
        vector<uint64_t> mu_i;
        if (i < proto_params.ell) {
            mu_i = m2_data[i];   // ct1 plaintext = 1, so slot i decodes to μ_2[i]
        } else {
            size_t r = i - proto_params.ell;
            mu_i.assign(N, 0);
            for (size_t c = 0; c < proto_params.ell; c++) {
                auto prod = mk::plain_mul(sk2.T_coeffs[r][c], m2_data[c], proto_params);
                mu_i = mk::plain_add(mu_i, prod, proto_params);
            }
        }
        // Build shift = (K−1)·μ_i·correction  mod p, then lift to R_q (NTT form).
        vector<uint64_t> shift_data(N);
        for (size_t j = 0; j < N; j++) {
            uint64_t v = mu_i[j];
            v = multiply_uint_mod(v, K_minus_1_mod_p, proto_params.p);
            v = multiply_uint_mod(v, tct.correction, proto_params.p);
            shift_data[j] = v;
        }
        mk::RnsPoly shift_ntt = mk::plain_to_rns_ntt(shift_data, proto_params, tct.d3[i].num_primes);
        tct.d3[i] = mk::add(tct.d3[i], shift_ntt, proto_params);

        // Add √K Gaussian noise to d3[i] only (no s-amplification in mkDec).
        mk::add_gaussian_p_noise(tct.d3[i], sigma_K, proto_params, prng_proto);
    }
    long t8 = take_us();
    cout << "  Step 8 (sim 2^20 adds → K·μ, √K·σ)    " << fmt(t8)
         << "       est ~" << fixed << setprecision(1)
         << std::log2(sigma_K * std::sqrt(2.0 * std::log((double)N))) << "\n";

    // Step 9 SKIPPED — non-ZK variant. No smudge applied to the tensor.
    cout << "  Step 9 (SKIPPED — no smudge)              0       —\n";

    // Step 10: mod-switch tensor L=3 → L=2 → L=1.  With intermediate mkDec check.
    {
        t_start = steady_clock::now();
        mk::mod_switch_tensor(tct, proto_params);  // L=3 → L=2 (drop q2=48b)
        long t10a = take_us();
        cout << "  Step 10a (mod-switch L=3 → L=2)       " << fmt(t10a) << "\n";

        // Intermediate mkDec at L=2 to isolate where it breaks.
        {
            bool sp;
            auto dec_L2 = mk::mkDec(sk1_proto, sk2, tct, proto_params, sp);
            bool match = true;
            for (size_t i = 0; i < proto_params.ell; i++) {
                vector<uint64_t> expected(N);
                uint64_t K_mp = SKIP_STEP_8 ? 1 : (K % p_val);
                for (size_t j = 0; j < N; j++)
                    expected[j] = multiply_uint_mod(m2_data[i][j], K_mp, proto_params.p);
                if (dec_L2[i] != expected) { match = false; break; }
            }
            cout << "  [intermediate mkDec @ L=2]            " << (match ? "OK" : "FAIL")
                 << " sparsi: " << (sp ? "ok" : "FAIL") << "\n";
        }

        t_start = steady_clock::now();
        mk::mod_switch_tensor(tct, proto_params);  // L=2 → L=1 (drop q1=49b)
        long t10b = take_us();
        cout << "  Step 10b (mod-switch L=2 → L=1)       " << fmt(t10b) << "\n";
    }

    cout << "\n";

    // ============================================================
    // Final mkDec & verify
    // ============================================================
    bool spchk;
    auto dec_out = mk::mkDec(sk1_proto, sk2, tct, proto_params, spchk);

    const uint64_t K_mod_p = SKIP_STEP_8 ? 1 : (K % p_val);
    bool all_match = true;
    size_t first_mismatch_slot = (size_t)-1, first_mismatch_idx = 0;
    uint64_t got_v = 0, want_v = 0;
    for (size_t i = 0; i < proto_params.ell; i++) {
        vector<uint64_t> expected(N);
        for (size_t j = 0; j < N; j++)
            expected[j] = multiply_uint_mod(m2_data[i][j], K_mod_p, proto_params.p);
        if (dec_out[i] != expected) {
            all_match = false;
            if (first_mismatch_slot == (size_t)-1) {
                first_mismatch_slot = i;
                for (size_t j = 0; j < N; j++) {
                    if (dec_out[i][j] != expected[j]) {
                        first_mismatch_idx = j;
                        got_v = dec_out[i][j];
                        want_v = expected[j];
                        break;
                    }
                }
            }
        }
    }
    cout << "[final mkDec on ℓ=" << proto_params.ell << " slots]     "
         << (all_match ? "OK" : "FAIL")
         << "  (expected: " << (SKIP_STEP_8 ? "μ_2[i]" : "2^20 · μ_2[i]") << " mod p)\n";
    if (!all_match) {
        cout << "  first mismatch: slot " << first_mismatch_slot
             << ", coeff " << first_mismatch_idx
             << ": got " << got_v << ", want " << want_v << "\n";
    }
    cout << "[sparsification check]            " << (spchk ? "OK" : "FAIL") << "\n";

    // Final ciphertext size: 12 R_q polys at q0 (33 bits) × 8192 coefficients each.
    const int q0_final_bits = mk::bit_length(proto_params.q_chain[0].value());
    const size_t tensor_polys = 1 + tct.d1.size() + 1 + tct.d3.size();
    const size_t bytes = tensor_polys * N * q0_final_bits / 8;
    cout << "\n[final ciphertext at L=1]"
         << "\n  components:   " << tensor_polys << " R_q polys (2 + 2·ell')"
         << "\n  coefficients: " << N << " per poly"
         << "\n  coeff width:  " << q0_final_bits << " bits"
         << "\n  total size:   " << bytes << " bytes (" << bytes / 1024.0 << " KiB)\n";

    cout << "\n[initial noise budget needed]   " << budget_init << " bits\n";

    return all_match ? 0 : 1;
}
