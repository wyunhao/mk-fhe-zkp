// SEAL standard BGV + BFV benchmark.
//   N = 8192, t (plain modulus) = 7340033 (= 7·2^20 + 1)
//   data-level q ≈ 130 bits (chain {44, 43, 43} + 60-bit special prime)
//
// Times each op 100 iterations and reports the average.

#include "seal/seal.h"
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace seal;
using namespace std;
using namespace std::chrono;

static double bench(const string& label, int iters, std::function<void()> f) {
    f();  // warm-up
    auto t0 = steady_clock::now();
    for (int i = 0; i < iters; i++) f();
    auto t1 = steady_clock::now();
    double avg_ms = duration_cast<nanoseconds>(t1 - t0).count() / double(iters) / 1.0e6;
    cout << "  " << left << setw(28) << label
         << "  avg " << fixed << setprecision(3) << setw(8) << avg_ms << " ms\n";
    return avg_ms;
}

static void run_scheme(scheme_type scheme, const string& scheme_name) {
    const size_t N = 8192;
    const uint64_t t_val = (uint64_t(7) << 20) + 1;   // 7340033

    EncryptionParameters parms(scheme);
    parms.set_poly_modulus_degree(N);
    // {44, 43, 43, 60} → first_context_data has 3 data primes (130-bit q); 60b is special.
    parms.set_coeff_modulus(CoeffModulus::Create(N, {44, 43, 43, 60}));
    parms.set_plain_modulus(t_val);

    SEALContext ctx(parms, /*expand_mod_chain*/ true, sec_level_type::none);

    auto& cd = *ctx.first_context_data();
    int total_bits = 0;
    for (auto& m : cd.parms().coeff_modulus()) {
        int b = 0; uint64_t v = m.value(); while (v) { b++; v >>= 1; }
        total_bits += b;
    }

    cout << "\n=========================================================\n";
    cout << "=== SEAL " << scheme_name << " benchmark ===\n";
    cout << "  N = " << N << ", t = " << t_val
         << ",  data-level log q = " << total_bits << " bits"
         << ",  sec_level = none\n";
    cout << "=========================================================\n";

    KeyGenerator kg(ctx);
    SecretKey sk = kg.secret_key();
    PublicKey pk; kg.create_public_key(pk);
    RelinKeys rk; kg.create_relin_keys(rk);
    GaloisKeys gk; kg.create_galois_keys(gk);

    Encryptor enc(ctx, pk);
    enc.set_secret_key(sk);
    Decryptor dec(ctx, sk);
    Evaluator eval(ctx);
    BatchEncoder batcher(ctx);

    vector<uint64_t> data(N);
    for (size_t i = 0; i < N; i++) data[i] = i % t_val;
    Plaintext pt_a; batcher.encode(data, pt_a);
    Plaintext pt_b = pt_a;

    Ciphertext ct_a_L3, ct_b_L3;
    enc.encrypt(pt_a, ct_a_L3);
    enc.encrypt(pt_b, ct_b_L3);

    // L=2 versions (one mod-switch dropped)
    Ciphertext ct_a_L2 = ct_a_L3; eval.mod_switch_to_next_inplace(ct_a_L2);
    Ciphertext ct_b_L2 = ct_b_L3; eval.mod_switch_to_next_inplace(ct_b_L2);

    const int iters = 100;

    cout << "\n-- L=3 (full chain, log q = " << total_bits << " bits) --\n";
    bench("EncryptSecret",        iters, [&]{ Ciphertext c; enc.encrypt_symmetric(pt_a, c); });
    bench("EncryptPublic",        iters, [&]{ Ciphertext c; enc.encrypt(pt_a, c); });
    bench("Decrypt",              iters, [&]{ Plaintext p; dec.decrypt(ct_a_L3, p); });
    bench("EncodeBatch",          iters, [&]{ Plaintext p; batcher.encode(data, p); });
    bench("DecodeBatch",          iters, [&]{ vector<uint64_t> d; batcher.decode(pt_a, d); });
    bench("EvaluateAddCt",        iters, [&]{ Ciphertext c; eval.add(ct_a_L3, ct_b_L3, c); });
    bench("EvaluateAddCtInplace", iters, [&]{ Ciphertext c = ct_a_L3; eval.add_inplace(c, ct_b_L3); });
    bench("EvaluateNegate",       iters, [&]{ Ciphertext c; eval.negate(ct_a_L3, c); });
    bench("EvaluateSubCt",        iters, [&]{ Ciphertext c; eval.sub(ct_a_L3, ct_b_L3, c); });
    bench("EvaluateAddPt",        iters, [&]{ Ciphertext c = ct_a_L3; eval.add_plain_inplace(c, pt_b); });
    bench("EvaluateMulPt",        iters, [&]{ Ciphertext c; eval.multiply_plain(ct_a_L3, pt_b, c); });
    bench("EvaluateMulCt",        iters, [&]{ Ciphertext c; eval.multiply(ct_a_L3, ct_b_L3, c); });
    bench("EvaluateMulCtInplace", iters, [&]{ Ciphertext c = ct_a_L3; eval.multiply_inplace(c, ct_b_L3); });
    bench("EvaluateSquare",       iters, [&]{ Ciphertext c; eval.square(ct_a_L3, c); });
    bench("EvaluateModSwitchInpl",iters, [&]{ Ciphertext c = ct_a_L3; eval.mod_switch_to_next_inplace(c); });
    bench("EvaluateRelinInplace", iters, [&]{
        Ciphertext c; eval.multiply(ct_a_L3, ct_b_L3, c);
        eval.relinearize_inplace(c, rk);
    });
    bench("EvaluateRotateRows",   iters, [&]{ Ciphertext c; eval.rotate_rows(ct_a_L3, 1, gk, c); });
    bench("EvaluateRotateCols",   iters, [&]{ Ciphertext c; eval.rotate_columns(ct_a_L3, gk, c); });

    cout << "\n-- L=2 (after one mod-switch) --\n";
    bench("EvaluateAddCt",        iters, [&]{ Ciphertext c; eval.add(ct_a_L2, ct_b_L2, c); });
    bench("EvaluateAddCtInplace", iters, [&]{ Ciphertext c = ct_a_L2; eval.add_inplace(c, ct_b_L2); });
    bench("EvaluateMulPt",        iters, [&]{ Ciphertext c; eval.multiply_plain(ct_a_L2, pt_b, c); });
    bench("EvaluateMulCt",        iters, [&]{ Ciphertext c; eval.multiply(ct_a_L2, ct_b_L2, c); });
    bench("EvaluateMulCtInplace", iters, [&]{ Ciphertext c = ct_a_L2; eval.multiply_inplace(c, ct_b_L2); });
    bench("EvaluateSquare",       iters, [&]{ Ciphertext c; eval.square(ct_a_L2, c); });
}

int main() {
    cout << "Microsoft SEAL version " << SEAL_VERSION
         << " — BGV + BFV benchmarks at matching prototype params\n";
    run_scheme(scheme_type::bgv, "BGV");
    run_scheme(scheme_type::bfv, "BFV");
    return 0;
}
