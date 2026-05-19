// SEAL standard-BGV benchmark at matching params (N=8192, p=7340033, q={44,43,43}).
// Times individual ops 100x and reports averages, comparable to prototype.cpp `bench`.
#include "seal/seal.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace seal;
using namespace std;
using namespace std::chrono;

static double bench(const string& label, int iters, std::function<void()> f) {
    // Warm-up
    f();
    auto t0 = steady_clock::now();
    for (int i = 0; i < iters; i++) f();
    auto t1 = steady_clock::now();
    double avg_ms = duration_cast<nanoseconds>(t1 - t0).count() / double(iters) / 1.0e6;
    cout << "  " << left << setw(28) << label
         << "  avg " << fixed << setprecision(3) << setw(8) << avg_ms << " ms"
         << "  (over " << iters << " iters)\n";
    return avg_ms;
}

int main() {
    const size_t N = 8192;
    const uint64_t p_val = (uint64_t(7) << 20) + 1;  // 7340033

    EncryptionParameters parms(scheme_type::bgv);
    parms.set_poly_modulus_degree(N);
    // Match prototype's data-level q chain {44, 43, 43}. SEAL reserves the LAST prime in the
    // user-supplied list as the "special" key-switching prime, so we pass {44,43,43,special}.
    // first_context_data then sees 3 primes = 130 bits, matching our prototype's L=3.
    parms.set_coeff_modulus(CoeffModulus::Create(N, {44, 43, 43, 60}));
    parms.set_plain_modulus(p_val);

    SEALContext ctx(parms, /*expand_mod_chain*/ true, sec_level_type::none);

    auto& cd = *ctx.first_context_data();
    int total_bits = 0;
    for (auto& m : cd.parms().coeff_modulus()) {
        int b = 0; uint64_t v = m.value(); while (v) { b++; v >>= 1; }
        total_bits += b;
    }
    size_t L_first = cd.parms().coeff_modulus().size();
    cout << "=== SEAL BGV benchmark @ matching params ===\n"
         << "  N = " << N << ", p = " << p_val
         << ", data-level primes = " << L_first << ", log q = " << total_bits << " bits\n";
    cout << "  (key level adds an extra special prime for key-switching)\n";
    cout << "  scheme = BGV, sec_level = none\n\n";

    KeyGenerator kg(ctx);
    SecretKey sk = kg.secret_key();
    PublicKey pk; kg.create_public_key(pk);
    RelinKeys rk; kg.create_relin_keys(rk);

    Encryptor enc(ctx, pk);
    enc.set_secret_key(sk);
    Decryptor dec(ctx, sk);
    Evaluator eval(ctx);
    BatchEncoder batcher(ctx);

    // Make some plaintext + ciphertexts at L=3 (full chain).
    vector<uint64_t> data(N);
    for (size_t i = 0; i < N; i++) data[i] = i % p_val;
    Plaintext pt; batcher.encode(data, pt);

    Ciphertext ct_a_L3, ct_b_L3;
    enc.encrypt(pt, ct_a_L3);
    enc.encrypt(pt, ct_b_L3);

    // L=2 copies via one mod-switch.
    Ciphertext ct_a_L2 = ct_a_L3; eval.mod_switch_to_next_inplace(ct_a_L2);
    Ciphertext ct_b_L2 = ct_b_L3; eval.mod_switch_to_next_inplace(ct_b_L2);

    const int iters = 100;

    cout << "-- L=3 (full chain, log q = " << total_bits << ") --\n";
    bench("EncryptSecret",     iters, [&]{ Ciphertext c; enc.encrypt_symmetric(pt, c); });
    bench("EncryptPublic",     iters, [&]{ Ciphertext c; enc.encrypt(pt, c); });
    bench("Decrypt",           iters, [&]{ Plaintext p; dec.decrypt(ct_a_L3, p); });
    bench("AddCt",             iters, [&]{ Ciphertext c; eval.add(ct_a_L3, ct_b_L3, c); });
    bench("AddCtInplace",      iters, [&]{ Ciphertext c = ct_a_L3; eval.add_inplace(c, ct_b_L3); });
    bench("MulCt",             iters, [&]{ Ciphertext c; eval.multiply(ct_a_L3, ct_b_L3, c); });
    bench("MulCtInplace",      iters, [&]{ Ciphertext c = ct_a_L3; eval.multiply_inplace(c, ct_b_L3); });
    bench("Square",            iters, [&]{ Ciphertext c; eval.square(ct_a_L3, c); });
    bench("ModSwitchInplace",  iters, [&]{ Ciphertext c = ct_a_L3; eval.mod_switch_to_next_inplace(c); });
    bench("RelinInplace",      iters, [&]{
        Ciphertext c; eval.multiply(ct_a_L3, ct_b_L3, c);
        eval.relinearize_inplace(c, rk);
    });

    int dropped_bits = 0;
    {
        uint64_t v = cd.parms().coeff_modulus().back().value();
        while (v) { dropped_bits++; v >>= 1; }
    }
    cout << "\n-- L=2 (after one mod-switch, log q = " << (total_bits - dropped_bits) << ") --\n";
    bench("AddCt",             iters, [&]{ Ciphertext c; eval.add(ct_a_L2, ct_b_L2, c); });
    bench("AddCtInplace",      iters, [&]{ Ciphertext c = ct_a_L2; eval.add_inplace(c, ct_b_L2); });
    bench("MulCt",             iters, [&]{ Ciphertext c; eval.multiply(ct_a_L2, ct_b_L2, c); });
    bench("MulCtInplace",      iters, [&]{ Ciphertext c = ct_a_L2; eval.multiply_inplace(c, ct_b_L2); });
    bench("Square",            iters, [&]{ Ciphertext c; eval.square(ct_a_L2, c); });

    return 0;
}
