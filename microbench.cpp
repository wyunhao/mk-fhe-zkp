// Pin down where our mkMul time goes:
//   (A) raw arithmetic floor: 24 dyadic_product_coeffmod calls into preallocated buffers (8 muls × 3 primes)
//   (B) same but allocating 8 fresh "RnsPoly"-shaped buffers per iteration (mirrors mkMul)
//   (C) SEAL's MulCt for reference
#include "seal/seal.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/ntt.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace seal;
using namespace seal::util;
using namespace std;
using namespace std::chrono;

static double bench(const string& label, int iters, std::function<void()> f) {
    f();  // warmup
    auto t0 = steady_clock::now();
    for (int i = 0; i < iters; i++) f();
    auto t1 = steady_clock::now();
    double avg_ms = duration_cast<nanoseconds>(t1 - t0).count() / double(iters) / 1.0e6;
    cout << "  " << left << setw(50) << label
         << "  avg " << fixed << setprecision(3) << setw(8) << avg_ms << " ms\n";
    return avg_ms;
}

int main() {
    const size_t N = 8192;
    const uint64_t p_val = (uint64_t(7) << 20) + 1;

    EncryptionParameters parms(scheme_type::bgv);
    parms.set_poly_modulus_degree(N);
    parms.set_coeff_modulus(CoeffModulus::Create(N, {44, 43, 43, 60}));  // last = special prime
    parms.set_plain_modulus(p_val);
    SEALContext ctx(parms, true, sec_level_type::none);
    auto& cd = *ctx.first_context_data();
    auto coeff_modulus = cd.parms().coeff_modulus();
    const size_t L = coeff_modulus.size();
    auto& kd = *ctx.key_context_data();
    cout << "  key_context_data primes: " << kd.parms().coeff_modulus().size()
         << ", first_context_data primes: " << L << "\n";

    cout << "=== Microbench: where does mkMul time go? ===\n";
    cout << "  N = " << N << ", L = " << L << " (full chain ~130 bits)\n\n";

    // Preallocated input buffers (NTT form, dummy data).
    vector<vector<uint64_t>> A(L, vector<uint64_t>(N));
    vector<vector<uint64_t>> B(L, vector<uint64_t>(N));
    for (size_t k = 0; k < L; k++) {
        uint64_t q = coeff_modulus[k].value();
        for (size_t i = 0; i < N; i++) {
            A[k][i] = (uint64_t(i) * 1000003ULL + k) % q;
            B[k][i] = (uint64_t(i) * 2000003ULL + k * 7) % q;
        }
    }

    const int iters = 200;

    // (A) Raw arithmetic floor: 8 dyadic products per prime, preallocated output bufs.
    vector<vector<uint64_t>> OUT8(8, vector<uint64_t>(L * N));
    bench("(A) 8 dyadic_products/prime, preallocated output ", iters, [&]{
        for (size_t k = 0; k < L; k++) {
            for (int m = 0; m < 8; m++) {
                dyadic_product_coeffmod(A[k].data(), B[k].data(), N, coeff_modulus[k], OUT8[m].data() + k*N);
            }
        }
    });

    // (B) Same arithmetic but allocating 8 fresh L*N-sized buffers each call (mimics mkMul).
    bench("(B) 8 muls/prime + 8 fresh L*N alloc+zero each call ", iters, [&]{
        for (int m = 0; m < 8; m++) {
            vector<uint64_t> buf(L * N, 0ULL);   // allocate AND zero — same as make_zero()
            for (size_t k = 0; k < L; k++) {
                dyadic_product_coeffmod(A[k].data(), B[k].data(), N, coeff_modulus[k], buf.data() + k*N);
            }
        }
    });

    // (B') Same as (B) but skip zeroing (since we overwrite anyway).
    bench("(B') 8 muls/prime + 8 fresh L*N alloc (no zero)     ", iters, [&]{
        for (int m = 0; m < 8; m++) {
            vector<uint64_t> buf; buf.resize(L * N);  // no zero
            for (size_t k = 0; k < L; k++) {
                dyadic_product_coeffmod(A[k].data(), B[k].data(), N, coeff_modulus[k], buf.data() + k*N);
            }
        }
    });

    // (C) SEAL's MulCt for reference (4 muls + 1 add per prime, tiled).
    KeyGenerator kg(ctx);
    PublicKey pk; kg.create_public_key(pk);
    Encryptor enc(ctx, pk);
    Evaluator eval(ctx);
    BatchEncoder batcher(ctx);
    vector<uint64_t> data(N);
    for (size_t i = 0; i < N; i++) data[i] = i % p_val;
    Plaintext pt; batcher.encode(data, pt);
    Ciphertext ca, cb;
    enc.encrypt(pt, ca);
    enc.encrypt(pt, cb);
    bench("(C) SEAL MulCt (4 muls + 1 add / prime, tiled)      ", iters, [&]{
        Ciphertext c; eval.multiply(ca, cb, c);
    });

    // (D) 4 dyadic_products + 1 add per prime, preallocated — what an "optimal mkMul-skeleton" would
    // hit if it did SEAL-equivalent work (and there were no multi-key fan-out).
    vector<vector<uint64_t>> OUT4(4, vector<uint64_t>(L * N));
    vector<vector<uint64_t>> TMP(L, vector<uint64_t>(N));
    bench("(D) 4 muls + 1 add per prime, preallocated          ", iters, [&]{
        for (size_t k = 0; k < L; k++) {
            dyadic_product_coeffmod(A[k].data(), B[k].data(), N, coeff_modulus[k], OUT4[0].data() + k*N);
            dyadic_product_coeffmod(A[k].data(), B[k].data(), N, coeff_modulus[k], TMP[k].data());
            dyadic_product_coeffmod(A[k].data(), B[k].data(), N, coeff_modulus[k], OUT4[1].data() + k*N);
            add_poly_coeffmod(OUT4[1].data() + k*N, TMP[k].data(), N, coeff_modulus[k], OUT4[1].data() + k*N);
            dyadic_product_coeffmod(A[k].data(), B[k].data(), N, coeff_modulus[k], OUT4[2].data() + k*N);
        }
    });

    return 0;
}
