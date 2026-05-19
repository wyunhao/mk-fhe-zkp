// Extended-field multi-key BGV mul: ct_2 encrypts (R_p[Y]/f(Y))^{ell'} with d = deg f(Y) = 5.
// ct_1 stays as a normal BGV ciphertext encrypting an element of R_p; its plaintext is treated
// as a constant in Y (only the Y^0 coefficient is nonzero).
//
// In this construction, ciphertext components live in (R_p[Y]/f(Y))_q ≅ (R_q)^d:
//   ct_1 = (a_1, b_1) with a_1, b_1 in R_q       (constant in Y)
//   ct_2 = (a(Y), b_0(Y), ..., b_{ell'-1}(Y))
//        where a(Y), b_i(Y) are elements of (R_p[Y]/f(Y))_q, each represented by d R_q-polys.
//
// Because ct_1's components are constants in Y, mkMul never invokes the f(Y) reduction:
//   d0(Y)   [k] = a_1 * a_2[k]          for k in [d]              (d products)
//   d2(Y)   [k] = b_1 * a_2[k]          for k in [d]              (d products)
//   d1[i](Y)[k] = a_1 * b_2[i][k]       for i in [ell'], k in [d] (ell'*d products)
//   d3[i](Y)[k] = b_1 * b_2[i][k]       for i in [ell'], k in [d] (ell'*d products)
// Total dyadic products per prime: 2d (shared a-side) + 2*ell'*d (per-slot) = 2d(1+ell').
// With ell' = 5 and d = 5: 2*5*6 = 60 dyadic products per prime (3 primes => 180 at L=3).
//
// This file does the BENCHMARK only. The mkMul kernel is the same operation our main prototype
// already validates end-to-end; lifting to (R_p[Y]/f(Y))^{ell'} multiplies the per-prime work by
// a factor of d (and adds the 2d shared a-side terms). No new arithmetic is introduced because
// the f(Y) reduction never fires when one operand is Y-constant.

#include <chrono>
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

using namespace seal;
using namespace seal::util;
using namespace std;
using namespace std::chrono;

namespace ext {

struct Params {
    size_t N = 0;
    int log_N = 0;
    Modulus p;
    vector<Modulus> q_chain;
    vector<unique_ptr<NTTTables>> q_ntt;
};

static void params_init(Params& params, uint64_t plain_val, const vector<int>& q_bits, size_t N) {
    params.N = N;
    params.log_N = 0;
    while ((size_t(1) << params.log_N) < N) params.log_N++;
    if ((size_t(1) << params.log_N) != N) throw runtime_error("N must be power of 2");

    params.p = Modulus(plain_val);
    for (int bs : q_bits) {
        auto cands = get_primes(2 * uint64_t(N), bs, 8);
        Modulus chosen(0);
        for (auto& c : cands) {
            bool dup = false;
            for (auto& q : params.q_chain) if (q.value() == c.value()) { dup = true; break; }
            if (!dup) { chosen = c; break; }
        }
        if (chosen.value() == 0) throw runtime_error("could not pick distinct q prime");
        params.q_chain.push_back(chosen);
    }
    for (auto& q : params.q_chain) {
        params.q_ntt.emplace_back(make_unique<NTTTables>(params.log_N, q));
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

static RnsPoly sample_uniform(const Params& params, size_t L, mt19937_64& rng) {
    RnsPoly r = make_zero(params, L);
    for (size_t k = 0; k < L; k++) {
        uint64_t q = params.q_chain[k].value();
        uint64_t* dst = row(r, k);
        for (size_t i = 0; i < params.N; i++) dst[i] = rng() % q;
    }
    return r;
}

static RnsPoly mul(const RnsPoly& a, const RnsPoly& b, const Params& params) {
    RnsPoly r = make_zero(params, a.num_primes);
    for (size_t k = 0; k < a.num_primes; k++) {
        dyadic_product_coeffmod(row(a, k), row(b, k), a.N, params.q_chain[k], row(r, k));
    }
    return r;
}

// ---- ciphertext types ----
// ct_1: standard BGV scalar ct in R_q (same as prototype.cpp).
struct CipherScalar { RnsPoly a, b; };

// ct_2: vector ct over the extension (R_p[Y]/f(Y))^{ell'}.
// Components live in (R_p[Y]/f(Y))_q, each represented by d R_q-polys (one per Y-coefficient).
struct CipherVectorExt {
    vector<RnsPoly> a;            // length d (one R_q-poly per Y-coefficient)
    vector<vector<RnsPoly>> b;    // [ell'][d]
};

// Tensor output of mkMul on (CipherScalar, CipherVectorExt). For each i in [ell'] and k in [d]
// we get four R_q tensor components (d0, d2, d1[i], d3[i])[k]. d0[k], d2[k] are shared across i.
struct TensorCtExt {
    vector<RnsPoly> d0;            // [d]
    vector<RnsPoly> d2;            // [d]
    vector<vector<RnsPoly>> d1;    // [ell'][d]
    vector<vector<RnsPoly>> d3;    // [ell'][d]
};

static TensorCtExt mkMul_ext(const CipherScalar& ct1, const CipherVectorExt& ct2, const Params& params) {
    const size_t d = ct2.a.size();
    const size_t ellp = ct2.b.size();
    TensorCtExt out;

    out.d0.reserve(d);
    out.d2.reserve(d);
    for (size_t k = 0; k < d; k++) {
        out.d0.emplace_back(mul(ct1.a, ct2.a[k], params));     // a_1 * a_2[k]
        out.d2.emplace_back(mul(ct1.b, ct2.a[k], params));     // b_1 * a_2[k]
    }

    out.d1.resize(ellp);
    out.d3.resize(ellp);
    for (size_t i = 0; i < ellp; i++) {
        out.d1[i].reserve(d);
        out.d3[i].reserve(d);
        for (size_t k = 0; k < d; k++) {
            out.d1[i].emplace_back(mul(ct1.a, ct2.b[i][k], params));   // a_1 * b_2[i][k]
            out.d3[i].emplace_back(mul(ct1.b, ct2.b[i][k], params));   // b_1 * b_2[i][k]
        }
    }
    return out;
}

// Drop the last prime from an RnsPoly (simulates the post-mod-switch level for benchmarking).
// We don't reconstruct correct decryption here — only the number/shape of dyadic products matters
// for runtime measurement.
static RnsPoly truncate_last(const RnsPoly& r, const Params& params) {
    RnsPoly out;
    out.N = params.N;
    out.num_primes = r.num_primes - 1;
    out.data.assign(r.data.begin(), r.data.begin() + (r.num_primes - 1) * params.N);
    (void)params;
    return out;
}

} // namespace ext

int main() {
    using namespace ext;

    Params params;
    const uint64_t p_val = (uint64_t(7) << 20) + 1;   // 7*2^20 + 1
    const size_t N = 8192;
    const size_t ell = 4, tau = 1, ell_prime = ell + tau;   // ell' = 5
    const size_t d = 5;                                      // deg f(Y)
    const vector<int> q_bits = {44, 43, 43};

    params_init(params, p_val, q_bits, N);

    int total_bits = 0;
    for (auto& q : params.q_chain) {
        int b = 0; uint64_t v = q.value(); while (v) { b++; v >>= 1; }
        total_bits += b;
    }
    int last_bits = 0;
    { uint64_t v = params.q_chain.back().value(); while (v) { last_bits++; v >>= 1; } }

    cout << "=== Extended multi-key BGV mul (ct_2 over (R_p[Y]/f(Y))^{ell'}) ===\n";
    cout << "  N = " << N << ", p = " << p_val << "\n";
    cout << "  ell = " << ell << ", tau = " << tau << ", ell' = " << ell_prime
         << ",  d = deg f(Y) = " << d << "\n";
    cout << "  log q @ L=3 = " << total_bits << " bits"
         << ",  log q @ L=2 = " << (total_bits - last_bits) << " bits\n";
    cout << "  dyadic products per prime in mkMul_ext = 2d(1 + ell') = "
         << 2 * d * (1 + ell_prime) << "\n";
    cout << "  (no f(Y) reduction in mkMul since ct_1 is constant in Y)\n\n";

    mt19937_64 rng(12345);

    // Build dummy ciphertexts at L=3. (Random uniform polys in NTT form. We only care about
    // dyadic-product timing; we're not verifying decryption here.)
    CipherScalar ct1_L3{ sample_uniform(params, 3, rng), sample_uniform(params, 3, rng) };

    CipherVectorExt ct2_L3;
    ct2_L3.a.reserve(d);
    for (size_t k = 0; k < d; k++) ct2_L3.a.emplace_back(sample_uniform(params, 3, rng));
    ct2_L3.b.resize(ell_prime);
    for (size_t i = 0; i < ell_prime; i++) {
        ct2_L3.b[i].reserve(d);
        for (size_t k = 0; k < d; k++) ct2_L3.b[i].emplace_back(sample_uniform(params, 3, rng));
    }

    // L=2 versions by truncating the last prime.
    CipherScalar ct1_L2{ truncate_last(ct1_L3.a, params), truncate_last(ct1_L3.b, params) };
    CipherVectorExt ct2_L2;
    ct2_L2.a.reserve(d);
    for (auto& a : ct2_L3.a) ct2_L2.a.emplace_back(truncate_last(a, params));
    ct2_L2.b.resize(ell_prime);
    for (size_t i = 0; i < ell_prime; i++) {
        ct2_L2.b[i].reserve(d);
        for (auto& b : ct2_L3.b[i]) ct2_L2.b[i].emplace_back(truncate_last(b, params));
    }

    auto bench = [](const string& label, int iters, function<void()> f) {
        f();  // warm-up
        auto t0 = steady_clock::now();
        for (int i = 0; i < iters; i++) f();
        auto t1 = steady_clock::now();
        double ms = duration_cast<nanoseconds>(t1 - t0).count() / double(iters) / 1.0e6;
        cout << "  " << left << setw(22) << label
             << "  avg " << fixed << setprecision(3) << setw(8) << ms << " ms"
             << "  (over " << iters << " iters)\n";
    };

    const int iters = 100;
    uint64_t sink = 0;
    bench("mkMul_ext @ L=3", iters, [&]{
        auto r = mkMul_ext(ct1_L3, ct2_L3, params);
        sink ^= r.d0[0].data[0] ^ r.d3.back().back().data[0];
    });
    bench("mkMul_ext @ L=2", iters, [&]{
        auto r = mkMul_ext(ct1_L2, ct2_L2, params);
        sink ^= r.d0[0].data[0] ^ r.d3.back().back().data[0];
    });
    (void)sink;

    return 0;
}
