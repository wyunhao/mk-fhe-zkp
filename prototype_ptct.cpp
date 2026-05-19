// Plaintext × Ciphertext "scalar multiplication" test for standard BGV (single-key).
//
// Two flavors covered:
//   (A) INTEGER scalar c × ct:   c ∈ [0, p), broadcast to all coefficients.
//                                Noise grows deterministically: ||e'||_inf = c · ||e||_inf.
//   (B) POLY  scalar m × ct:     m ∈ R_p (a plaintext polynomial in NTT form).
//                                Noise grows as ||m · e||_inf ≤ N · ||m||_inf · ||e||_inf
//                                (worst case) or sqrt(N·log N) · ||m||_canon · σ_e (typical).
//
// Both run on a fresh ct at L=2 (q_chain = {44, 43} → 87-bit q) so noise inf-norm can be
// reconstructed exactly via uint128 CRT — same machinery as in the main prototype.
//
// Reports:
//   - log2(||e||_inf) before/after each variant
//   - runtime per multiplication (averaged over 1000 iters)
//   - end-to-end decryption check against (c·μ) or (m·μ) mod p

#include <chrono>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "seal/modulus.h"
#include "seal/util/ntt.h"
#include "seal/util/numth.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/uintarithsmallmod.h"

using namespace seal;
using namespace seal::util;
using namespace std;
using namespace std::chrono;

namespace ptct {

// -------------------------- Params, RnsPoly --------------------------
struct Params {
    size_t N = 0;
    int log_N = 0;
    Modulus p;
    vector<Modulus> q_chain;
    vector<unique_ptr<NTTTables>> q_ntt;
    int sigma_eta = 21;
    size_t hw = 32;
};

static int bit_length(uint64_t v) { int b = 0; while (v) { b++; v >>= 1; } return b; }

static void params_init(Params& p, uint64_t plain_val, const vector<int>& q_bits, size_t N) {
    p.N = N;
    p.log_N = 0;
    while ((size_t(1) << p.log_N) < N) p.log_N++;
    if ((size_t(1) << p.log_N) != N) throw runtime_error("N must be power of 2");
    p.p = Modulus(plain_val);
    for (int bs : q_bits) {
        auto cands = get_primes(2 * uint64_t(N), bs, 8);
        Modulus chosen(0);
        for (auto& c : cands) {
            bool dup = false;
            for (auto& q : p.q_chain) if (q.value() == c.value()) { dup = true; break; }
            if (!dup) { chosen = c; break; }
        }
        if (chosen.value() == 0) throw runtime_error("could not pick distinct q prime");
        p.q_chain.push_back(chosen);
    }
    for (auto& q : p.q_chain) p.q_ntt.emplace_back(make_unique<NTTTables>(p.log_N, q));
}

struct RnsPoly {
    vector<uint64_t> data;
    size_t num_primes = 0;
    size_t N = 0;
};

static inline uint64_t* row(RnsPoly& r, size_t k) { return r.data.data() + k * r.N; }
static inline const uint64_t* row(const RnsPoly& r, size_t k) { return r.data.data() + k * r.N; }

static RnsPoly make_zero(const Params& p, size_t L) {
    RnsPoly r; r.num_primes = L; r.N = p.N; r.data.assign(L * p.N, 0ULL); return r;
}

static RnsPoly sample_uniform(const Params& p, size_t L, mt19937_64& rng) {
    RnsPoly r = make_zero(p, L);
    for (size_t k = 0; k < L; k++) {
        uint64_t q = p.q_chain[k].value();
        uint64_t* dst = row(r, k);
        for (size_t i = 0; i < p.N; i++) dst[i] = rng() % q;
    }
    return r;
}

static RnsPoly small_signed_to_ntt(const vector<int>& coeffs, const Params& p, size_t L) {
    RnsPoly r = make_zero(p, L);
    for (size_t k = 0; k < L; k++) {
        uint64_t q = p.q_chain[k].value();
        uint64_t* dst = row(r, k);
        for (size_t i = 0; i < p.N; i++) {
            int v = coeffs[i];
            dst[i] = (v >= 0) ? uint64_t(v) : (q - uint64_t(-v));
        }
        ntt_negacyclic_harvey(dst, *p.q_ntt[k]);
    }
    return r;
}

static RnsPoly plain_to_rns_ntt(const vector<uint64_t>& coeffs, const Params& p, size_t L) {
    RnsPoly r = make_zero(p, L);
    for (size_t k = 0; k < L; k++) {
        memcpy(row(r, k), coeffs.data(), p.N * sizeof(uint64_t));
        ntt_negacyclic_harvey(row(r, k), *p.q_ntt[k]);
    }
    return r;
}

static int cbd(int eta, mt19937_64& rng) {
    int s = 0;
    for (int i = 0; i < eta; i++) s += int(rng() & 1ULL);
    for (int i = 0; i < eta; i++) s -= int(rng() & 1ULL);
    return s;
}

static RnsPoly add(const RnsPoly& a, const RnsPoly& b, const Params& p) {
    RnsPoly r = make_zero(p, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++)
        add_poly_coeffmod(row(a, k), row(b, k), a.N, p.q_chain[k], row(r, k));
    return r;
}
static RnsPoly sub(const RnsPoly& a, const RnsPoly& b, const Params& p) {
    RnsPoly r = make_zero(p, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++)
        sub_poly_coeffmod(row(a, k), row(b, k), a.N, p.q_chain[k], row(r, k));
    return r;
}
static RnsPoly mul_dyadic(const RnsPoly& a, const RnsPoly& b, const Params& p) {
    RnsPoly r = make_zero(p, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++)
        dyadic_product_coeffmod(row(a, k), row(b, k), a.N, p.q_chain[k], row(r, k));
    return r;
}
static RnsPoly mul_by_p(const RnsPoly& a, const Params& p) {
    RnsPoly r = make_zero(p, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++) {
        uint64_t p_mod_qk = barrett_reduce_64(p.p.value(), p.q_chain[k]);
        multiply_poly_scalar_coeffmod(row(a, k), a.N, p_mod_qk, p.q_chain[k], row(r, k));
    }
    return r;
}

// -------------------------- BGV scheme (single-key) --------------------------
struct CipherScalar { RnsPoly a, b; };

static RnsPoly keygen(const Params& p, mt19937_64& rng) {
    vector<int> coeffs(p.N, 0);
    vector<size_t> positions(p.N);
    for (size_t i = 0; i < p.N; i++) positions[i] = i;
    for (size_t i = 0; i < p.hw; i++) {
        size_t j = i + (rng() % (p.N - i));
        swap(positions[i], positions[j]);
    }
    for (size_t i = 0; i < p.hw; i++) coeffs[positions[i]] = (rng() & 1ULL) ? 1 : -1;
    return small_signed_to_ntt(coeffs, p, p.q_chain.size());
}

static CipherScalar encrypt(const RnsPoly& s_ntt, const vector<uint64_t>& msg,
                            const Params& p, mt19937_64& rng) {
    const size_t L = p.q_chain.size();
    CipherScalar ct;
    ct.a = sample_uniform(p, L, rng);
    vector<int> e_signed(p.N);
    for (size_t i = 0; i < p.N; i++) e_signed[i] = cbd(p.sigma_eta, rng);
    RnsPoly e_ntt = small_signed_to_ntt(e_signed, p, L);
    RnsPoly as = mul_dyadic(ct.a, s_ntt, p);
    RnsPoly pe = mul_by_p(e_ntt, p);
    RnsPoly mu_ntt = plain_to_rns_ntt(msg, p, L);
    ct.b = add(add(as, pe, p), mu_ntt, p);
    return ct;
}

// L=2 noise + decryption (uint128 CRT reconstruction).
// Given expected_mu in [0, p), returns (decoded msg in [0, p), log2 of noise inf-norm).
struct DecResult { vector<uint64_t> msg; double log2_noise_inf; bool consistent; };

static DecResult decrypt_with_noise(const RnsPoly& s_ntt, const CipherScalar& ct,
                                    const vector<uint64_t>& expected_mu, const Params& p) {
    if (ct.a.num_primes != 2) throw runtime_error("decrypt_with_noise: requires L=2");
    RnsPoly as = mul_dyadic(ct.a, s_ntt, p);
    RnsPoly bma = sub(ct.b, as, p);
    vector<uint64_t> c0(p.N), c1(p.N);
    memcpy(c0.data(), row(bma, 0), p.N * sizeof(uint64_t));
    memcpy(c1.data(), row(bma, 1), p.N * sizeof(uint64_t));
    inverse_ntt_negacyclic_harvey(c0.data(), *p.q_ntt[0]);
    inverse_ntt_negacyclic_harvey(c1.data(), *p.q_ntt[1]);

    const uint64_t q0 = p.q_chain[0].value();
    const uint64_t q1 = p.q_chain[1].value();
    const __uint128_t Q = __uint128_t(q0) * __uint128_t(q1);
    const __uint128_t Q_half = Q >> 1;
    uint64_t inv_q0_mod_q1;
    try_invert_uint_mod(q0, p.q_chain[1], inv_q0_mod_q1);
    const int64_t p_val = (int64_t)p.p.value();

    DecResult res;
    res.msg.resize(p.N);
    res.consistent = true;
    __int128_t max_abs = 0;

    for (size_t i = 0; i < p.N; i++) {
        uint64_t r0 = c0[i], r1 = c1[i];
        uint64_t r0_mod_q1 = barrett_reduce_64(r0, p.q_chain[1]);
        uint64_t y = (r1 >= r0_mod_q1) ? (r1 - r0_mod_q1) : (q1 - (r0_mod_q1 - r1));
        y = multiply_uint_mod(y, inv_q0_mod_q1, p.q_chain[1]);
        __uint128_t x = __uint128_t(r0) + __uint128_t(q0) * __uint128_t(y);
        __int128_t v = (x > Q_half) ? (__int128_t)x - (__int128_t)Q : (__int128_t)x;

        int64_t mu_e = (int64_t)expected_mu[i];
        __int128_t diff = v - mu_e;
        if (diff % p_val != 0) {
            __int128_t diff2 = v - (mu_e - p_val);
            if (diff2 % p_val != 0) { res.consistent = false; }
            else diff = diff2;
        }
        __int128_t e = diff / p_val;
        if (e < 0) e = -e;
        if (e > max_abs) max_abs = e;

        // decode m: centered v mod p
        int64_t r = (int64_t)(v % p_val);
        if (r < 0) r += p_val;
        res.msg[i] = (uint64_t)r;
    }
    if (max_abs == 0) res.log2_noise_inf = 0.0;
    else {
        int hi = 0;
        __uint128_t u = (__uint128_t)max_abs;
        while (u > 1) { u >>= 1; hi++; }
        res.log2_noise_inf = hi + log2((double)(__uint128_t)max_abs / ((__uint128_t)1 << hi));
    }
    return res;
}

// -------------------------- Scalar multiplications --------------------------
// (A) Integer scalar: c × ct. c < p < q_0.
static void multiply_int_inplace(CipherScalar& ct, uint64_t c, const Params& p) {
    for (size_t k = 0; k < ct.a.num_primes; k++) {
        uint64_t c_qk = barrett_reduce_64(c, p.q_chain[k]);
        multiply_poly_scalar_coeffmod(row(ct.a, k), ct.a.N, c_qk, p.q_chain[k], row(ct.a, k));
        multiply_poly_scalar_coeffmod(row(ct.b, k), ct.b.N, c_qk, p.q_chain[k], row(ct.b, k));
    }
}

// (B) Polynomial scalar: m × ct.  m is given as an RNS-NTT poly already.
static CipherScalar multiply_plain(const CipherScalar& ct, const RnsPoly& m_ntt, const Params& p) {
    CipherScalar out;
    out.a = mul_dyadic(ct.a, m_ntt, p);
    out.b = mul_dyadic(ct.b, m_ntt, p);
    return out;
}

// -------------------------- Plaintext-side arithmetic --------------------------
// (used to compute the *expected* result μ·m mod p)
static unique_ptr<NTTTables> make_p_ntt(const Params& p) {
    return make_unique<NTTTables>(p.log_N, p.p);
}

static vector<uint64_t> plain_mul_mod_p(const vector<uint64_t>& a, const vector<uint64_t>& b,
                                        const Params& p, const NTTTables& p_ntt) {
    vector<uint64_t> an(a), bn(b), c(p.N);
    ntt_negacyclic_harvey(an.data(), p_ntt);
    ntt_negacyclic_harvey(bn.data(), p_ntt);
    dyadic_product_coeffmod(an.data(), bn.data(), p.N, p.p, c.data());
    inverse_ntt_negacyclic_harvey(c.data(), p_ntt);
    return c;
}

} // namespace ptct

// ============================================================================
int main() {
    using namespace ptct;

    Params p;
    const uint64_t p_val = (uint64_t(7) << 20) + 1;   // 7340033
    params_init(p, p_val, /*q_bits*/ {44, 43}, /*N*/ 8192);
    auto p_ntt = make_p_ntt(p);

    cout << "=== BGV plaintext × ciphertext scalar mul test ===\n"
         << "  N = " << p.N << ", p = " << p_val
         << ",  q chain {" << bit_length(p.q_chain[0].value()) << "b, "
                            << bit_length(p.q_chain[1].value()) << "b}"
         << " = " << (bit_length(p.q_chain[0].value()) + bit_length(p.q_chain[1].value())) << "b\n"
         << "  ternary HW = " << p.hw << ", CBD eta = " << p.sigma_eta << "\n\n";

    mt19937_64 rng(20260519);
    auto sk = keygen(p, rng);

    // Random plaintext for ct
    vector<uint64_t> mu(p.N);
    for (size_t i = 0; i < p.N; i++) mu[i] = rng() % p_val;

    auto ct_fresh = encrypt(sk, mu, p, rng);
    auto dec_fresh = decrypt_with_noise(sk, ct_fresh, mu, p);
    cout << fixed << setprecision(3);
    cout << "[fresh ct]                          log2||e||_inf = " << dec_fresh.log2_noise_inf
         << "   decrypt: " << (dec_fresh.msg == mu ? "OK" : "FAIL") << "\n\n";

    // ---------------- (A) Integer scalar c × ct ----------------
    {
        uint64_t c = (rng() % (p_val - 1)) + 1;   // c ∈ [1, p)
        CipherScalar ct = ct_fresh;
        multiply_int_inplace(ct, c, p);

        // expected μ_out = (c * mu) mod p
        vector<uint64_t> mu_out(p.N);
        for (size_t i = 0; i < p.N; i++) mu_out[i] = multiply_uint_mod(mu[i], c, p.p);

        auto dec_a = decrypt_with_noise(sk, ct, mu_out, p);
        cout << "(A) Integer scalar c × ct  (c = " << c
             << ", log2 c = " << log2((double)c) << "b):\n"
             << "    pre-mul   log2||e||_inf = " << dec_fresh.log2_noise_inf << "\n"
             << "    post-mul  log2||e||_inf = " << dec_a.log2_noise_inf
             << "   (predicted: pre + log2(c) = "
             << (dec_fresh.log2_noise_inf + log2((double)c)) << "b)\n"
             << "    decrypt: " << (dec_a.msg == mu_out ? "OK" : "FAIL") << "\n";

        // Time it
        const int iters = 1000;
        CipherScalar warm = ct_fresh;
        multiply_int_inplace(warm, c, p);
        auto t0 = steady_clock::now();
        for (int i = 0; i < iters; i++) {
            CipherScalar w = ct_fresh;
            multiply_int_inplace(w, c, p);
            asm volatile("" : : "r"(w.b.data[0]));  // prevent DCE
        }
        auto t1 = steady_clock::now();
        double us = duration_cast<nanoseconds>(t1 - t0).count() / double(iters) / 1.0e3;
        cout << "    runtime:  " << us << " us / op  (over " << iters << " iters)\n\n";
    }

    // ---------------- (B) Polynomial scalar m × ct ----------------
    {
        vector<uint64_t> m_coeffs(p.N);
        for (size_t i = 0; i < p.N; i++) m_coeffs[i] = rng() % p_val;
        RnsPoly m_ntt = plain_to_rns_ntt(m_coeffs, p, /*L=*/2);

        auto ct_out = multiply_plain(ct_fresh, m_ntt, p);
        auto mu_out = plain_mul_mod_p(mu, m_coeffs, p, *p_ntt);
        auto dec_b = decrypt_with_noise(sk, ct_out, mu_out, p);

        // Predicted noise (canonical-norm extension): sqrt(N · log N) · ||m||_canon · σ_e
        //   ||m||_canon  ≈ sqrt(N) · σ_m_unif = sqrt(N) · p/sqrt(12)
        //   σ_e          ≈ sqrt(η/2) ≈ 3.24  (fresh CBD)
        // ⇒  predicted ≈ sqrt(N) · (p/sqrt(12)) · σ_e · sqrt(log N) · O(1)
        double log2_pred = 0.5 * log2((double)p.N) + log2((double)p_val) - 0.5 * log2(12.0)
                          + log2(sqrt((double)p.sigma_eta / 2.0)) + 0.5 * log2(log((double)p.N));
        cout << "(B) Polynomial scalar m × ct  (m: random uniform in R_p):\n"
             << "    pre-mul   log2||e||_inf = " << dec_fresh.log2_noise_inf << "\n"
             << "    post-mul  log2||e||_inf = " << dec_b.log2_noise_inf
             << "   (canonical-norm prediction ≈ " << log2_pred << "b)\n"
             << "    decrypt: " << (dec_b.msg == mu_out ? "OK" : "FAIL") << "\n";

        const int iters = 1000;
        auto warm = multiply_plain(ct_fresh, m_ntt, p);
        auto t0 = steady_clock::now();
        for (int i = 0; i < iters; i++) {
            auto w = multiply_plain(ct_fresh, m_ntt, p);
            asm volatile("" : : "r"(w.b.data[0]));
        }
        auto t1 = steady_clock::now();
        double us = duration_cast<nanoseconds>(t1 - t0).count() / double(iters) / 1.0e3;
        cout << "    runtime:  " << us << " us / op  (over " << iters
             << " iters; excludes NTT-encoding the plaintext m)\n";
    }

    return 0;
}
