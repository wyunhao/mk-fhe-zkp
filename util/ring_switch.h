// Ring-switching utility: convert a BGV ciphertext at N_big = 16384 to N_small = 8192.
//
// Strategy: the secret key sk_big in R_{N_big} is constructed as the embedding of an underlying
// sk_small ∈ R_{N_small} via the map X → X² (so sk_big[2i] = sk_small[i], sk_big[2i+1] = 0).
// This makes sk_big fixed by the Galois automorphism σ: X → X^{N_big+1} ≡ -X (mod X^{N_big}+1),
// because σ fixes X² and hence every polynomial in X² is invariant.
//
// The trace map Tr(p) := p + σ(p) maps any polynomial in R_{N_big} to one whose odd-degree
// coefficients are zero (because σ negates X and the odd-degree terms cancel) and whose even-
// degree coefficients are doubled. Applied to the ciphertext components and to the encoded
// plaintext, the trace lifts the decryption identity from R_{N_big} into the σ-fixed subring,
// which is canonically isomorphic to R_{N_small} via X² ↔ Y.
//
// Procedure:
//   1. ct_trace = ct + apply_galois(ct, k = N_big + 1)
//   2. Extract every-other-coefficient from ct_trace.a, ct_trace.b → a_small, b_small ∈ R_{N_small}.
//   3. Decrypt under sk_small: (b_small + a_small · sk_small) mod q — yields 2μ + p · noise
//      (the trace doubled the encoded plaintext).
//   4. Multiply by 2⁻¹ mod p to recover μ.
//
// Why this works as a "key-switch" between sk_big and sk_small without a separate key-switching
// key: we deliberately *constructed* sk_big as the embedding of sk_small, so the Galois key for
// σ (which SEAL provides natively) doubles as the ring-switching key. If one wanted an arbitrary
// (independently sampled) sk_big, a true BV/HPS key-switching key from sk_big to embed(sk_small)
// would be needed before applying the trace.

#pragma once

#include <cstdint>
#include <vector>

#include "seal/seal.h"

namespace ring_switch_util {

// -------------------------- Keys --------------------------

// Sample a ternary secret with fixed Hamming weight in R_N (returns N signed coefficients).
std::vector<int> sample_ternary_hw(std::size_t N, std::size_t hw, std::uint64_t seed);

// Embed sk_small ∈ R_{N_small} into R_{2·N_small} via X → X²:
//   sk_big[2i] = sk_small[i],  sk_big[2i+1] = 0.
std::vector<int> embed_sk(const std::vector<int>& sk_small);

// Build a SEAL SecretKey from a coefficient-form polynomial (length = poly_modulus_degree).
// Internally writes the coefficients into the key context's prime slabs and NTT-transforms.
seal::SecretKey make_seal_secret_key(const std::vector<int>& sk_coeffs,
                                     const seal::SEALContext& ctx);

// Create the Galois key needed for ring switching (k = N_big + 1, the X → -X automorphism).
seal::GaloisKeys make_ring_switch_galois_key(seal::KeyGenerator& kg, std::size_t N_big);

// -------------------------- BV-style key-switching key with gadget decomposition --------------
// Switches a ciphertext from sk_src to sk_dst (both polynomials in R_{N_big}, at the same RNS
// chain). Implemented in BV style: RNS decomposition (one branch per data prime) *combined*
// with a gadget decomposition base B = 2^log_base inside each branch.
//
// For each (k, b) with k ∈ [0, L) (RNS prime) and b ∈ [0, num_digits), the KSK has one piece
// encrypting sk_src · T_k · B^b under sk_dst. Apply decomposes each coefficient of (c_1 mod q_k)
// into base-B digits and accumulates digit · piece (k, b).
//
// Noise per piece is ≈ B · σ_e · √(N log N), so total post-apply noise ≈ √(L · num_digits) · B
// · σ_e · √(N log N).  For N=16384, L=3, B=256 (log_base=8), q_k=55 b → num_digits=7:
//   ≈ √21 · 256 · 3.24 · √(2·16384·14) ≈ 2^20 inf-norm  (≈ +20 b of noise).
// Smaller B → smaller per-piece noise but more pieces (more storage and more apply work).
struct KSwitchKey {
    struct Piece {
        std::vector<std::vector<std::uint64_t>> a;   // [L][N_big], NTT form per prime
        std::vector<std::vector<std::uint64_t>> b;   // [L][N_big], NTT form per prime
    };
    std::vector<std::vector<Piece>> pieces;          // pieces[k][b], k=RNS prime, b=base-B digit
    std::vector<seal::Modulus> q_chain;
    std::size_t N_big = 0;
    std::uint64_t plain_modulus = 0;
    int log_base = 8;                                // gadget base = 2^log_base
    int num_digits = 0;                              // pieces per RNS prime
};

// Build KSK from sk_src to sk_dst at the level corresponding to `parms_id`. Both sk's are given
// as signed coefficient vectors of length N_big. log_base sets the gadget base B=2^log_base.
KSwitchKey gen_key_switch_key(const std::vector<int>& sk_src_coeffs,
                              const std::vector<int>& sk_dst_coeffs,
                              const seal::SEALContext& ctx,
                              const seal::parms_id_type& parms_id,
                              std::uint64_t plain_modulus,
                              int log_base,
                              std::uint64_t seed);

// Apply the KSK to a SEAL Ciphertext (assumed under sk_src) — produces a fresh SEAL Ciphertext
// under sk_dst at the same level. Internally constructs the new ct in NTT form by RNS
// decomposition of c_1 (the "a" component in SEAL convention is c_1).
seal::Ciphertext apply_key_switch(const seal::Ciphertext& ct,
                                  const KSwitchKey& ksk,
                                  const seal::SEALContext& ctx);

// -------------------------- Plaintext helpers --------------------------

// Pack data_small (N_small values in [0, p)) into raw R_{N_big} coefficients at even indices;
// returns the N_big-length vector ready for a Plaintext.
std::vector<std::uint64_t> embed_plaintext(const std::vector<std::uint64_t>& data_small,
                                           std::size_t N_big);

// Build a SEAL Plaintext from raw coefficients.
seal::Plaintext make_seal_plaintext(const std::vector<std::uint64_t>& coeffs);

// -------------------------- Ring switching --------------------------

// A ciphertext in R_{N_small} at the same RNS chain as the input. Stored coefficient-wise
// (NOT in NTT form) per prime so decryption is straightforward.
struct RingSwitchedCt {
    // a[k][i] = i-th coefficient of "a polynomial" mod q_k.
    std::vector<std::vector<std::uint64_t>> a;
    std::vector<std::vector<std::uint64_t>> b;
    std::vector<seal::Modulus> q_chain;   // RNS primes at the ciphertext's current level
    std::size_t N_small = 0;
    std::uint64_t correction = 1;         // SEAL's BGV correction_factor at this level
    int log_q_total = 0;                  // for reporting
};

// Apply ring switching to a SEAL Ciphertext.  Assumes sk_big is the embedding of an sk_small
// (so the Galois automorphism with element N_big+1 is in `gk_ringswitch`).
RingSwitchedCt ring_switch(const seal::Ciphertext& ct_big,
                           const seal::Evaluator& eval,
                           const seal::GaloisKeys& gk_ringswitch,
                           const seal::SEALContext& ctx_big,
                           std::size_t N_small);

// Decrypt a RingSwitchedCt using sk_small.  Returns the recovered plaintext as N_small values
// in [0, p).  Internally compensates for (a) the doubling done by the trace map and (b) any
// mod-switching correction.
std::vector<std::uint64_t> decrypt_ring_switched(const RingSwitchedCt& rs_ct,
                                                 const std::vector<int>& sk_small_coeffs,
                                                 std::uint64_t plain_modulus);

// Decrypt and ALSO report the noise inf-norm (in log2 bits) of the ring-switched ciphertext.
// Same as decrypt_ring_switched but additionally writes the noise estimate into *out_log2_noise.
// Self-measures noise from the decrypted plaintext — no `expected_data` needed.
std::vector<std::uint64_t> decrypt_and_measure_noise(const RingSwitchedCt& rs_ct,
                                                     const std::vector<int>& sk_small_coeffs,
                                                     std::uint64_t plain_modulus,
                                                     double* out_log2_noise);

}  // namespace ring_switch_util
