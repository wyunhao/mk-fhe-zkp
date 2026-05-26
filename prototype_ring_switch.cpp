// Ring-switching test driver.
//
// Pipeline:
//   1. Initialize a BGV ciphertext at N=16384 with ~250 b initial noise budget.
//      (top-level data q chain {55, 55, 55, 55, 55} + 60 special = 335 b → 250 b budget.)
//   2. Mod-switch to reduce log q below 218 b, so once we're at N=8192 we're 128-bit-secure.
//   3. Ring-switch (trace + extract) the resulting ciphertext from R_{16384} to R_{8192}.
//   4. Decrypt the result under the 8192-dimension secret key and compare with the original
//      8192 raw plaintext slots.

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "seal/seal.h"
#include "util/ring_switch.h"

using namespace seal;
using namespace std;
using namespace std::chrono;

static int bit_length(uint64_t v) { int b = 0; while (v) { b++; v >>= 1; } return b; }

int main() {
    const size_t N_big   = 16384;
    const size_t N_small = 8192;
    const uint64_t p_val = (uint64_t(7) << 20) + 1;   // 7340033

    // ============================================================
    // Phase 1: SEAL setup at N=16384 with ~250 b initial budget.
    // ============================================================
    EncryptionParameters parms(scheme_type::bgv);
    parms.set_poly_modulus_degree(N_big);
    parms.set_coeff_modulus(CoeffModulus::Create(N_big, {55, 55, 55, 55, 55, 60}));
    parms.set_plain_modulus(p_val);
    // sec_level_type::none because we'll mod-switch *down* into the 128-bit-secure regime; SEAL's
    // initial check rejects the top chain since it's chosen specifically to give us headroom
    // before mod-switching. The post-mod-switch + ring-switched final ct *is* 128-bit secure.
    SEALContext ctx(parms, /*expand_mod_chain*/ true, sec_level_type::none);

    cout << "=== Ring-switch test: N=16384 → N=8192 ===\n";
    cout << "  SEAL chain at N_big=" << N_big << " (data primes + special):  ";
    int top_bits = 0;
    {
        auto& key_mods = ctx.key_context_data()->parms().coeff_modulus();
        for (size_t i = 0; i + 1 < key_mods.size(); i++) {
            int b = bit_length(key_mods[i].value());
            cout << b << "b ";
            top_bits += b;
        }
        cout << "+ " << bit_length(key_mods.back().value()) << "b(special)";
    }
    cout << "\n  Top-level data log q = " << top_bits << "b (target ~250 b init budget)\n";

    // ============================================================
    // Phase 2: TWO independent secret keys, then a REAL key-switching key from sk_big → sk_target.
    //   sk_big      ← random ternary HW=32 in R_{N_big}.
    //   sk_small    ← random ternary HW=32 in R_{N_small}.
    //   sk_target   = embed(sk_small), i.e., sk_small lifted to R_{N_big} via X→X² (invariant
    //                 under the X→-X Galois auto so the trace works).
    //   KSK         = BV-style RNS key-switching key from sk_big → sk_target.
    // The Galois key for the X→-X automorphism is built under sk_target (so it can be applied
    // *after* the key-switch puts the ct under sk_target).
    // ============================================================
    const uint64_t sk_big_seed   = 11111;
    const uint64_t sk_small_seed = 22222;
    const uint64_t ksk_seed      = 33333;

    auto sk_big_coeffs    = ring_switch_util::sample_ternary_hw(N_big,   /*hw*/ 32, sk_big_seed);
    auto sk_small_coeffs  = ring_switch_util::sample_ternary_hw(N_small, /*hw*/ 32, sk_small_seed);
    auto sk_target_coeffs = ring_switch_util::embed_sk(sk_small_coeffs);

    SecretKey sk_big_seal    = ring_switch_util::make_seal_secret_key(sk_big_coeffs,    ctx);
    SecretKey sk_target_seal = ring_switch_util::make_seal_secret_key(sk_target_coeffs, ctx);

    KeyGenerator kg_big(ctx, sk_big_seal);
    PublicKey pk;
    kg_big.create_public_key(pk);

    // Galois key under sk_target — used AFTER the key-switch so the ct is under sk_target.
    KeyGenerator kg_target(ctx, sk_target_seal);
    auto gk_ringswitch = ring_switch_util::make_ring_switch_galois_key(kg_target, N_big);

    Encryptor enc(ctx, pk);
    enc.set_secret_key(sk_big_seal);
    Decryptor dec_big(ctx, sk_big_seal);
    Decryptor dec_target(ctx, sk_target_seal);     // for noise inspection after the KSK
    Evaluator eval(ctx);

    // ============================================================
    // Phase 3: Random raw plaintext — 8192 values in [0, p), embedded at even positions in R_big.
    // ============================================================
    mt19937_64 rng(424242);
    vector<uint64_t> data_small(N_small);
    for (size_t i = 0; i < N_small; i++) data_small[i] = rng() % p_val;
    auto pt_coeffs = ring_switch_util::embed_plaintext(data_small, N_big);
    Plaintext pt = ring_switch_util::make_seal_plaintext(pt_coeffs);

    Ciphertext ct;
    enc.encrypt_symmetric(pt, ct);
    cout << "\n[step 1] Encrypted at top level.\n";
    cout << "  Initial noise budget for ct: " << dec_big.invariant_noise_budget(ct) << " b\n";

    // Sanity: SEAL decrypt should recover the embedded polynomial.
    {
        Plaintext pt_back;
        dec_big.decrypt(ct, pt_back);
        bool ok = (pt_back.coeff_count() <= N_big);
        for (size_t i = 0; i < N_big && ok; i++) {
            uint64_t got = (i < pt_back.coeff_count()) ? pt_back[i] : 0ULL;
            if (got != pt_coeffs[i]) { ok = false; break; }
        }
        cout << "  Sanity SEAL-decrypt at top: " << (ok ? "OK" : "FAIL") << "\n";
    }

    // ============================================================
    // Phase 4: Mod-switch to reduce log q below 218 b (so N_small=8192 is 128-bit-secure).
    // ============================================================
    cout << "\n[step 2] Mod-switching down: target log q < 218 b (128-bit secure at N=8192)\n"
            "         but keep ≥ 3 primes for now — KSK adds ~q_k·sqrt(N·logN)·σ_e per piece\n"
            "         which we'll mod-switch away AFTER applying the KSK.\n";
    while (ct.coeff_modulus_size() > 3) {
        auto cd_ptr = ctx.get_context_data(ct.parms_id());
        int cur_bits = 0;
        for (auto& m : cd_ptr->parms().coeff_modulus()) cur_bits += bit_length(m.value());
        cout << "  log q = " << cur_bits << "b (" << ct.coeff_modulus_size()
             << " primes) → mod-switch.\n";
        eval.mod_switch_to_next_inplace(ct);
    }
    {
        auto cd_ptr = ctx.get_context_data(ct.parms_id());
        int cur_bits = 0;
        for (auto& m : cd_ptr->parms().coeff_modulus()) cur_bits += bit_length(m.value());
        cout << "  Settled at log q = " << cur_bits << "b ("
             << ct.coeff_modulus_size() << " primes). Under 218 b → 128-bit-secure at N=" << N_small << ".\n";
    }
    cout << "  Budget at this level: " << dec_big.invariant_noise_budget(ct) << " b\n";

    // Verify SEAL ct still decrypts correctly.
    {
        Plaintext pt_back;
        dec_big.decrypt(ct, pt_back);
        bool ok = (pt_back.coeff_count() <= N_big);
        for (size_t i = 0; i < N_big && ok; i++) {
            uint64_t got = (i < pt_back.coeff_count()) ? pt_back[i] : 0ULL;
            if (got != pt_coeffs[i]) { ok = false; break; }
        }
        cout << "  Post-mod-switch SEAL-decrypt: " << (ok ? "OK" : "FAIL") << "\n";
    }

    // ============================================================
    // Phase 5: Build the KSK and run the FULL pipeline:
    //   (a) apply_key_switch(ct under sk_big)  → ct' under sk_target. Real noise from KSK.
    //   (b) ct'_sigma = apply_galois(ct', k=N_big+1) under sk_target.
    //   (c) ct_trace = ct' + ct'_sigma.
    //   (d) extract every-other-coefficient + decrypt under sk_small.
    // Budget reported at each stage.
    // ============================================================
    // Gadget base = 2^log_base.  Smaller log_base → less noise but more storage / apply work.
    //   log_base = 60 → no gadget (RNS-only, ~62 b noise).
    //   log_base = 16 → ~30 b noise, ~4 digits per prime.
    //   log_base = 8  → ~20 b noise, ~7 digits per prime  ← good default.
    //   log_base = 4  → ~17 b noise, ~14 digits per prime.
    int log_base = 8;
    if (const char* e = getenv("LOG_BASE")) log_base = atoi(e);

    cout << "\n[step 3] Building KSK (sk_big → sk_target = embed(sk_small)).\n";
    cout << "  gadget log_base = " << log_base << "  (base B = 2^" << log_base
         << " = " << (1ULL << log_base) << ")\n";
    auto t0 = steady_clock::now();
    auto ksk = ring_switch_util::gen_key_switch_key(sk_big_coeffs, sk_target_coeffs, ctx,
                                                     ct.parms_id(), p_val, log_base, ksk_seed);
    auto t1 = steady_clock::now();
    cout << "  gen_key_switch_key runtime: " << duration_cast<microseconds>(t1 - t0).count()
         << " us  (KSK has " << ksk.pieces.size() << " RNS primes × "
         << ksk.num_digits << " gadget digits = "
         << (ksk.pieces.size() * ksk.num_digits) << " pieces total)\n";

    int budget_pre = dec_big.invariant_noise_budget(ct);
    cout << "\n[step 4] Apply pipeline.\n";
    cout << "  budget before pipeline (under sk_big):              " << budget_pre << " b\n";

    // (a) Key-switch from sk_big to sk_target.  This is the big noise hit (BV-style KSK).
    auto t_ksw0 = steady_clock::now();
    Ciphertext ct_target = ring_switch_util::apply_key_switch(ct, ksk, ctx);
    auto t_ksw1 = steady_clock::now();
    int budget_after_ksw = dec_target.invariant_noise_budget(ct_target);
    cout << "  after KSK (sk_big → sk_target):                     "
         << budget_after_ksw << " b   (Δ = " << (budget_after_ksw - budget_pre) << " b, "
         << duration_cast<microseconds>(t_ksw1 - t_ksw0).count() << " us)\n";

    // Sanity: SEAL-decrypt ct_target under sk_target. If the polynomial matches pt_coeffs,
    // the KSK is structurally correct (just noise-heavy).
    {
        Plaintext pt_back;
        try {
            dec_target.decrypt(ct_target, pt_back);
            bool ok = true;
            for (size_t i = 0; i < N_big && ok; i++) {
                uint64_t got = (i < pt_back.coeff_count()) ? pt_back[i] : 0ULL;
                if (got != pt_coeffs[i]) { ok = false; break; }
            }
            cout << "  [sanity] SEAL-decrypt of ct_target under sk_target: "
                 << (ok ? "OK (KSK structurally correct)" : "FAIL (KSK has a logic bug)") << "\n";
        } catch (const std::exception& e) {
            cout << "  [sanity] decrypt threw: " << e.what() << "\n";
        }
    }

    // (a+) Mod-switch DOWN once to scrub the KSK noise (drops it by ~q_dropped bits).
    //     Brings ct from 3 primes → 2 primes (so our CRT decrypt at L=2 works downstream).
    auto t_ksw_ms0 = steady_clock::now();
    eval.mod_switch_to_next_inplace(ct_target);
    auto t_ksw_ms1 = steady_clock::now();
    int budget_after_ksw_ms = dec_target.invariant_noise_budget(ct_target);
    cout << "  + mod-switch to scrub KSK noise (3 → 2 primes):    "
         << budget_after_ksw_ms << " b   (Δ = "
         << (budget_after_ksw_ms - budget_after_ksw) << " b, "
         << duration_cast<microseconds>(t_ksw_ms1 - t_ksw_ms0).count() << " us)\n";

    // (b) Apply Galois X → -X under sk_target.
    Ciphertext ct_sigma;
    auto t_gal0 = steady_clock::now();
    eval.apply_galois(ct_target, static_cast<uint32_t>(N_big + 1), gk_ringswitch, ct_sigma);
    auto t_gal1 = steady_clock::now();
    int budget_after_galois = dec_target.invariant_noise_budget(ct_sigma);
    cout << "  after apply_galois (under sk_target):               "
         << budget_after_galois << " b   (Δ = "
         << (budget_after_galois - budget_after_ksw_ms) << " b, "
         << duration_cast<microseconds>(t_gal1 - t_gal0).count() << " us)\n";

    // (c) ct_trace = ct_target + ct_sigma.
    Ciphertext ct_trace;
    auto t_add0 = steady_clock::now();
    eval.add(ct_target, ct_sigma, ct_trace);
    auto t_add1 = steady_clock::now();
    int budget_after_trace = dec_target.invariant_noise_budget(ct_trace);
    cout << "  after add (trace):                                  "
         << budget_after_trace << " b   (Δ = "
         << (budget_after_trace - budget_after_galois) << " b, "
         << duration_cast<microseconds>(t_add1 - t_add0).count() << " us)\n";

    // (d) Extract every-other-coefficient → R_{N_small}.  Reuse ring_switch (Galois + add +
    // extract); pass ct_target (the post-mod-switch ct under sk_target) and gk_ringswitch.
    auto t_ex0 = steady_clock::now();
    auto rs_ct = ring_switch_util::ring_switch(ct_target, eval, gk_ringswitch, ctx, N_small);
    auto t_ex1 = steady_clock::now();
    cout << "  ring_switch (galois + add + extract):               -- "
         << duration_cast<microseconds>(t_ex1 - t_ex0).count() << " us\n";
    cout << "  result log q:          " << rs_ct.log_q_total << "b (" << rs_ct.q_chain.size()
         << " primes)\n";
    cout << "  result N_small:        " << rs_ct.N_small << "\n";
    cout << "  result correction:     " << rs_ct.correction << "\n";

    // ============================================================
    // Phase 6: Decrypt + self-measure post-ring-switch noise inf-norm (no expected_data).
    // ============================================================
    cout << "\n[step 5] Decrypting under sk_small + self-measuring noise.\n";
    double log2_noise_after = -1.0;
    auto t_dec0 = steady_clock::now();
    auto decoded = ring_switch_util::decrypt_and_measure_noise(rs_ct, sk_small_coeffs, p_val,
                                                                &log2_noise_after);
    auto t_dec1 = steady_clock::now();
    cout << "  decrypt runtime:       " << duration_cast<microseconds>(t_dec1 - t_dec0).count() << " us\n";

    // Convert log2(noise_inf) to a budget-style number for direct comparison.
    // budget_eq = log2(q) − log2(p) − log2(||e||) − 1.
    double budget_after = (double)rs_ct.log_q_total - std::log2((double)p_val) - log2_noise_after - 1.0;
    cout << fixed << setprecision(2);
    cout << "  post-ring-switch log2(||e||_inf) ≈ " << log2_noise_after
         << "   (equivalent budget ≈ " << budget_after << " b)\n";
    cout << "  Δ budget across full ring switch: " << (budget_after - budget_pre) << " b\n";

    bool ok = (decoded.size() == N_small);
    size_t first_bad = N_small;
    for (size_t i = 0; i < N_small && ok; i++) {
        if (decoded[i] != data_small[i]) {
            ok = false; first_bad = i;
            cout << "  Mismatch at idx " << i << ": got " << decoded[i]
                 << ", want " << data_small[i] << "\n";
            break;
        }
    }
    cout << "  Final verification:    " << (ok ? "OK" : "FAIL") << "\n";

    cout << "\n[summary]\n";
    cout << "  Started at N=" << N_big << ", ~" << top_bits << "b log q (~"
         << (top_bits - 28) << "b init budget)\n";
    cout << "  Mod-switched to log q = " << rs_ct.log_q_total << "b (128-bit-secure at N="
         << N_small << ")\n";
    cout << "  Ring-switched to N=" << N_small << " under sk_small (32-HW ternary)\n";
    cout << "  Final ciphertext: 2 R_q polys × " << N_small << " × " << rs_ct.log_q_total
         << " bits = " << (2 * N_small * rs_ct.log_q_total / 8) << " bytes\n";

    return ok ? 0 : 1;
}
