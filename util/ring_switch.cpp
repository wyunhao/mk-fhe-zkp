// Implementation of util/ring_switch.h — see header for algorithm description.

#include "ring_switch.h"

#include <cstring>
#include <random>
#include <stdexcept>

#include "seal/memorymanager.h"
#include "seal/modulus.h"
#include "seal/util/ntt.h"
#include "seal/util/numth.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/uintarithsmallmod.h"

namespace ring_switch_util {

using namespace seal;
using namespace seal::util;

// -------------------------- Keys --------------------------

std::vector<int> sample_ternary_hw(std::size_t N, std::size_t hw, std::uint64_t seed) {
    if (hw > N) throw std::runtime_error("hw > N");
    std::mt19937_64 rng(seed);
    std::vector<int> coeffs(N, 0);
    std::vector<std::size_t> positions(N);
    for (std::size_t i = 0; i < N; i++) positions[i] = i;
    for (std::size_t i = 0; i < hw; i++) {
        std::size_t j = i + static_cast<std::size_t>(rng() % (N - i));
        std::swap(positions[i], positions[j]);
    }
    for (std::size_t i = 0; i < hw; i++) {
        coeffs[positions[i]] = (rng() & 1ULL) ? 1 : -1;
    }
    return coeffs;
}

std::vector<int> embed_sk(const std::vector<int>& sk_small) {
    std::vector<int> sk_big(2 * sk_small.size(), 0);
    for (std::size_t i = 0; i < sk_small.size(); i++) sk_big[2 * i] = sk_small[i];
    return sk_big;
}

SecretKey make_seal_secret_key(const std::vector<int>& sk_coeffs,
                               const SEALContext& ctx) {
    auto& key_cd = *ctx.key_context_data();
    auto& parms = key_cd.parms();
    auto& coeff_modulus = parms.coeff_modulus();
    const std::size_t coeff_count = parms.poly_modulus_degree();
    const std::size_t coeff_modulus_size = coeff_modulus.size();
    if (sk_coeffs.size() != coeff_count) {
        throw std::runtime_error("sk coefficient count mismatch");
    }

    SecretKey sk;
    sk.data().resize(util::mul_safe(coeff_count, coeff_modulus_size));
    sk.data().parms_id() = parms_id_zero;
    std::uint64_t* sk_ptr = sk.data().data();

    // Write coefficients per prime, signed → mod q_k.
    for (std::size_t k = 0; k < coeff_modulus_size; k++) {
        const std::uint64_t q = coeff_modulus[k].value();
        for (std::size_t i = 0; i < coeff_count; i++) {
            const int v = sk_coeffs[i];
            sk_ptr[k * coeff_count + i] = (v >= 0)
                ? static_cast<std::uint64_t>(v)
                : (q - static_cast<std::uint64_t>(-v));
        }
    }

    // NTT in place per prime against key-context NTT tables.
    const auto& ntt_tables = key_cd.small_ntt_tables();
    for (std::size_t k = 0; k < coeff_modulus_size; k++) {
        ntt_negacyclic_harvey(sk_ptr + k * coeff_count, ntt_tables[k]);
    }
    sk.data().parms_id() = key_cd.parms_id();
    return sk;
}

GaloisKeys make_ring_switch_galois_key(KeyGenerator& kg, std::size_t N_big) {
    // X → -X corresponds to galois_elt = N_big + 1 (since X^{N+1} ≡ -X mod X^N+1).
    const std::uint32_t galois_elt = static_cast<std::uint32_t>(N_big + 1);
    GaloisKeys gk;
    kg.create_galois_keys(std::vector<std::uint32_t>{galois_elt}, gk);
    return gk;
}

// -------------------------- BV-style RNS key-switching --------------------------

// Embed a signed coefficient vector into NTT-form per RNS prime at the given level.
// Returns sk_ntt[k][i] for k ∈ [0, L), i ∈ [0, N).
static std::vector<std::vector<std::uint64_t>> sk_signed_to_ntt(
    const std::vector<int>& sk,
    const std::vector<Modulus>& q_chain,
    const std::vector<NTTTables*>& ntt_tables,
    std::size_t N) {
    const std::size_t L = q_chain.size();
    std::vector<std::vector<std::uint64_t>> out(L, std::vector<std::uint64_t>(N, 0));
    for (std::size_t k = 0; k < L; k++) {
        const std::uint64_t q = q_chain[k].value();
        for (std::size_t i = 0; i < N; i++) {
            const int v = sk[i];
            out[k][i] = (v >= 0) ? static_cast<std::uint64_t>(v)
                                 : (q - static_cast<std::uint64_t>(-v));
        }
        ntt_negacyclic_harvey(out[k].data(), *ntt_tables[k]);
    }
    return out;
}

// CBD(η=21) noise sample.
static int cbd_sample(std::mt19937_64& rng, int eta = 21) {
    int s = 0;
    for (int i = 0; i < eta; i++) s += static_cast<int>(rng() & 1ULL);
    for (int i = 0; i < eta; i++) s -= static_cast<int>(rng() & 1ULL);
    return s;
}

KSwitchKey gen_key_switch_key(const std::vector<int>& sk_src_coeffs,
                              const std::vector<int>& sk_dst_coeffs,
                              const SEALContext& ctx,
                              const parms_id_type& parms_id,
                              std::uint64_t plain_modulus,
                              int log_base,
                              std::uint64_t seed) {
    auto cd_ptr = ctx.get_context_data(parms_id);
    if (!cd_ptr) throw std::runtime_error("gen_key_switch_key: bad parms_id");
    auto& cd = *cd_ptr;
    const auto& q_chain = cd.parms().coeff_modulus();
    const std::size_t L = q_chain.size();
    const std::size_t N = cd.parms().poly_modulus_degree();
    if (sk_src_coeffs.size() != N || sk_dst_coeffs.size() != N) {
        throw std::runtime_error("gen_key_switch_key: sk coefficient count != N");
    }
    if (log_base <= 0 || log_base > 60) log_base = 8;

    const auto& ntt_tables_ref = cd.small_ntt_tables();
    std::vector<NTTTables*> ntt_ptrs(L, nullptr);
    for (std::size_t k = 0; k < L; k++) {
        ntt_ptrs[k] = const_cast<NTTTables*>(&ntt_tables_ref[k]);
    }

    auto sk_src_ntt = sk_signed_to_ntt(sk_src_coeffs, q_chain, ntt_ptrs, N);
    auto sk_dst_ntt = sk_signed_to_ntt(sk_dst_coeffs, q_chain, ntt_ptrs, N);

    Modulus p_mod(plain_modulus);

    // Determine num_digits: ceil(log2(max q_k) / log_base).
    int max_log_q = 0;
    for (auto& q : q_chain) {
        int b = 0; std::uint64_t v = q.value(); while (v) { b++; v >>= 1; }
        if (b > max_log_q) max_log_q = b;
    }
    const int num_digits = (max_log_q + log_base - 1) / log_base;

    KSwitchKey ksk;
    ksk.q_chain.assign(q_chain.begin(), q_chain.end());
    ksk.N_big = N;
    ksk.plain_modulus = plain_modulus;
    ksk.log_base = log_base;
    ksk.num_digits = num_digits;
    ksk.pieces.assign(L, std::vector<KSwitchKey::Piece>(num_digits));

    std::mt19937_64 rng(seed);

    std::vector<std::uint64_t> p_mod_qj(L);
    for (std::size_t j = 0; j < L; j++) p_mod_qj[j] = barrett_reduce_64(plain_modulus, q_chain[j]);

    // Build pieces[k][b]: encrypts (sk_src · T_k · 2^{b·log_base}) under sk_dst.
    // In NTT slab-j form:
    //   - slab j == k: ((sk_src · 2^{b·log_base}) mod q_k)  ← scalar multiply sk_src_ntt[k]
    //   - slab j != k: 0  (because T_k is 0 mod q_j for j ≠ k)
    for (std::size_t k = 0; k < L; k++) {
        const std::uint64_t qk = q_chain[k].value();
        for (int b = 0; b < num_digits; b++) {
            auto& piece = ksk.pieces[k][b];
            piece.a.assign(L, std::vector<std::uint64_t>(N, 0));
            piece.b.assign(L, std::vector<std::uint64_t>(N, 0));

            // Single small CBD polynomial sampled once, then reduced mod each q_j.
            std::vector<int> e_kb_signed(N);
            for (std::size_t i = 0; i < N; i++) e_kb_signed[i] = cbd_sample(rng);

            // pow_2_kb = 2^{b · log_base} mod qk (used only at slab j == k).
            // For b · log_base < log2(qk) (typical case), it's < qk. Otherwise compute mod qk.
            std::uint64_t pow2_b = 1;
            for (int bit = 0; bit < b * log_base; bit++) {
                pow2_b = (pow2_b << 1);
                if (pow2_b >= qk) pow2_b -= qk;
            }

            for (std::size_t j = 0; j < L; j++) {
                const std::uint64_t qj = q_chain[j].value();

                // a[j]: uniform random mod q_j (NTT form ≡ coeff form for uniform).
                for (std::size_t i = 0; i < N; i++) piece.a[j][i] = rng() % qj;

                // e_kb mod q_j → NTT.
                std::vector<std::uint64_t> e_kb_ntt(N);
                for (std::size_t i = 0; i < N; i++) {
                    int v = e_kb_signed[i];
                    e_kb_ntt[i] = (v >= 0) ? static_cast<std::uint64_t>(v)
                                            : (qj - static_cast<std::uint64_t>(-v));
                }
                ntt_negacyclic_harvey(e_kb_ntt.data(), *ntt_ptrs[j]);

                // a · sk_dst.
                std::vector<std::uint64_t> a_sk_dst(N);
                dyadic_product_coeffmod(piece.a[j].data(), sk_dst_ntt[j].data(), N, q_chain[j],
                                        a_sk_dst.data());

                // pe = p · e_kb_ntt.
                std::vector<std::uint64_t> pe(N);
                multiply_poly_scalar_coeffmod(e_kb_ntt.data(), N, p_mod_qj[j], q_chain[j], pe.data());

                // b[j] = pe - a · sk_dst.
                sub_poly_coeffmod(pe.data(), a_sk_dst.data(), N, q_chain[j], piece.b[j].data());

                // Add (sk_src · 2^{b·log_base}) at slab k only.
                if (j == k) {
                    std::vector<std::uint64_t> sk_scaled(N);
                    multiply_poly_scalar_coeffmod(sk_src_ntt[k].data(), N, pow2_b, q_chain[k],
                                                  sk_scaled.data());
                    add_poly_coeffmod(piece.b[j].data(), sk_scaled.data(), N, q_chain[j],
                                      piece.b[j].data());
                }
            }
        }
    }

    return ksk;
}

Ciphertext apply_key_switch(const Ciphertext& ct,
                            const KSwitchKey& ksk,
                            const SEALContext& ctx) {
    auto cd_ptr = ctx.get_context_data(ct.parms_id());
    if (!cd_ptr) throw std::runtime_error("apply_key_switch: bad ct parms_id");
    auto& cd = *cd_ptr;
    const auto& q_chain = cd.parms().coeff_modulus();
    const std::size_t L = q_chain.size();
    const std::size_t N = cd.parms().poly_modulus_degree();
    if (!ct.is_ntt_form()) throw std::runtime_error("apply_key_switch: expected NTT-form ct");
    if (ct.size() != 2) throw std::runtime_error("apply_key_switch: expected size-2 ct");
    if (L != ksk.q_chain.size()) throw std::runtime_error("apply_key_switch: KSK chain mismatch");
    if (N != ksk.N_big) throw std::runtime_error("apply_key_switch: KSK N mismatch");

    const auto& ntt_tables_ref = cd.small_ntt_tables();
    std::vector<const NTTTables*> ntt_ptrs(L, nullptr);
    for (std::size_t k = 0; k < L; k++) ntt_ptrs[k] = &ntt_tables_ref[k];

    // c0[k][.] and c1[k][.] (NTT form per prime).
    auto load_slab = [&](const std::uint64_t* base, std::size_t k) {
        std::vector<std::uint64_t> v(N);
        std::memcpy(v.data(), base + k * N, N * sizeof(std::uint64_t));
        return v;
    };

    // Initialize new_c0 = c0, new_c1 = 0.
    std::vector<std::vector<std::uint64_t>> new_c0(L), new_c1(L);
    for (std::size_t k = 0; k < L; k++) {
        new_c0[k] = load_slab(ct.data(0), k);
        new_c1[k].assign(N, 0);
    }

    // Gadget decomposition base B = 2^log_base.  Each c_1 mod q_k coefficient is decomposed
    // into num_digits base-B digits, each in [0, B).  Apply: for (kdig, b), accumulate
    //   digit · KSK.pieces[kdig][b]  into new_c0 / new_c1.
    const int log_base = ksk.log_base;
    const int num_digits = ksk.num_digits;
    const std::uint64_t base_mask = (std::uint64_t(1) << log_base) - 1ULL;

    for (std::size_t kdig = 0; kdig < L; kdig++) {
        std::vector<std::uint64_t> c1_kdig_coeffs(N);
        std::memcpy(c1_kdig_coeffs.data(), ct.data(1) + kdig * N, N * sizeof(std::uint64_t));
        inverse_ntt_negacyclic_harvey(c1_kdig_coeffs.data(), *ntt_ptrs[kdig]);

        for (int b = 0; b < num_digits; b++) {
            // Extract base-B digit b of each coefficient.
            // (b · log_base might exceed 63 for very large q_k, but our chain has q_k ≤ 60 b so
            //  shifts stay within uint64 range.)
            std::vector<std::uint64_t> c1_kb_coeffs(N);
            const int shift = b * log_base;
            for (std::size_t i = 0; i < N; i++) {
                c1_kb_coeffs[i] = (shift < 64) ? ((c1_kdig_coeffs[i] >> shift) & base_mask) : 0ULL;
            }

            for (std::size_t j = 0; j < L; j++) {
                // c1_kb_coeffs[i] ∈ [0, B) — already smaller than q_j (B = 2^log_base, log_base ≤ 30).
                std::vector<std::uint64_t> c1_kb_modqj(c1_kb_coeffs);
                ntt_negacyclic_harvey(c1_kb_modqj.data(), *ntt_ptrs[j]);

                std::vector<std::uint64_t> tmp(N);
                dyadic_product_coeffmod(c1_kb_modqj.data(), ksk.pieces[kdig][b].a[j].data(), N,
                                        q_chain[j], tmp.data());
                add_poly_coeffmod(new_c1[j].data(), tmp.data(), N, q_chain[j], new_c1[j].data());

                dyadic_product_coeffmod(c1_kb_modqj.data(), ksk.pieces[kdig][b].b[j].data(), N,
                                        q_chain[j], tmp.data());
                add_poly_coeffmod(new_c0[j].data(), tmp.data(), N, q_chain[j], new_c0[j].data());
            }
        }
    }

    // Build the output Ciphertext.
    Ciphertext out(ctx, ct.parms_id());
    out.resize(ctx, ct.parms_id(), 2);
    out.is_ntt_form() = true;
    out.scale() = ct.scale();
    out.correction_factor() = ct.correction_factor();

    for (std::size_t k = 0; k < L; k++) {
        std::memcpy(out.data(0) + k * N, new_c0[k].data(), N * sizeof(std::uint64_t));
        std::memcpy(out.data(1) + k * N, new_c1[k].data(), N * sizeof(std::uint64_t));
    }
    return out;
}

// -------------------------- Plaintext helpers --------------------------

std::vector<std::uint64_t> embed_plaintext(const std::vector<std::uint64_t>& data_small,
                                           std::size_t N_big) {
    if (data_small.size() * 2 != N_big) {
        throw std::runtime_error("embed_plaintext: data_small.size() must equal N_big/2");
    }
    std::vector<std::uint64_t> coeffs(N_big, 0);
    for (std::size_t i = 0; i < data_small.size(); i++) coeffs[2 * i] = data_small[i];
    return coeffs;
}

Plaintext make_seal_plaintext(const std::vector<std::uint64_t>& coeffs) {
    Plaintext pt;
    pt.resize(coeffs.size());
    for (std::size_t i = 0; i < coeffs.size(); i++) pt[i] = coeffs[i];
    return pt;
}

// -------------------------- Ring switching --------------------------

RingSwitchedCt ring_switch(const Ciphertext& ct_big,
                           const Evaluator& eval,
                           const GaloisKeys& gk_ringswitch,
                           const SEALContext& ctx_big,
                           std::size_t N_small) {
    // 1. Trace: ct_trace = ct + σ(ct) where σ(X) = X^{N_big+1} = -X.
    auto cd_ptr = ctx_big.get_context_data(ct_big.parms_id());
    if (!cd_ptr) throw std::runtime_error("ring_switch: bad ct parms_id");
    auto& cd = *cd_ptr;
    const std::size_t N_big = cd.parms().poly_modulus_degree();
    if (N_big != 2 * N_small) {
        throw std::runtime_error("ring_switch: N_big must equal 2 * N_small");
    }
    const std::uint32_t galois_elt = static_cast<std::uint32_t>(N_big + 1);

    Ciphertext ct_sigma;
    eval.apply_galois(ct_big, galois_elt, gk_ringswitch, ct_sigma);

    Ciphertext ct_trace;
    eval.add(ct_big, ct_sigma, ct_trace);

    if (ct_trace.size() != 2) {
        throw std::runtime_error("ring_switch: expected size-2 ciphertext (relin needed if larger)");
    }
    if (!ct_trace.is_ntt_form()) {
        throw std::runtime_error("ring_switch: expected NTT-form (BGV default)");
    }

    // 2. Inverse-NTT each prime slab, extract every-other-coefficient, store in coeff form.
    auto& coeff_modulus = cd.parms().coeff_modulus();
    const std::size_t L = coeff_modulus.size();
    const auto& ntt_big = cd.small_ntt_tables();

    RingSwitchedCt out;
    out.N_small = N_small;
    out.q_chain.reserve(L);
    out.a.assign(L, std::vector<std::uint64_t>(N_small, 0));
    out.b.assign(L, std::vector<std::uint64_t>(N_small, 0));

    int log_q_total = 0;
    for (std::size_t k = 0; k < L; k++) {
        out.q_chain.push_back(coeff_modulus[k]);
        std::uint64_t v = coeff_modulus[k].value();
        while (v) { log_q_total++; v >>= 1; }

        std::vector<std::uint64_t> a_big(N_big), b_big(N_big);
        // SEAL stores c0 = b, c1 = a (decryption: c0 + c1·s).  Keep that convention here.
        std::memcpy(b_big.data(), ct_trace.data(0) + k * N_big, N_big * sizeof(std::uint64_t));
        std::memcpy(a_big.data(), ct_trace.data(1) + k * N_big, N_big * sizeof(std::uint64_t));

        inverse_ntt_negacyclic_harvey(a_big.data(), ntt_big[k]);
        inverse_ntt_negacyclic_harvey(b_big.data(), ntt_big[k]);

        for (std::size_t i = 0; i < N_small; i++) {
            out.a[k][i] = a_big[2 * i];
            out.b[k][i] = b_big[2 * i];
        }
    }
    out.log_q_total = log_q_total;
    out.correction = ct_big.correction_factor();
    return out;
}

std::vector<std::uint64_t> decrypt_ring_switched(const RingSwitchedCt& rs,
                                                 const std::vector<int>& sk_small_coeffs,
                                                 std::uint64_t plain_modulus) {
    const std::size_t N = rs.N_small;
    const std::size_t L = rs.q_chain.size();
    if (sk_small_coeffs.size() != N) {
        throw std::runtime_error("decrypt_ring_switched: sk_small size mismatch");
    }
    if (L != 1 && L != 2) {
        throw std::runtime_error("decrypt_ring_switched: only L in {1, 2} supported");
    }

    // Build NTT tables for N_small at each q_k.
    int log_N = 0; while ((std::size_t(1) << log_N) < N) log_N++;
    std::vector<std::unique_ptr<NTTTables>> ntt_small;
    ntt_small.reserve(L);
    for (std::size_t k = 0; k < L; k++) {
        ntt_small.emplace_back(std::make_unique<NTTTables>(log_N, rs.q_chain[k],
                                                            MemoryManager::GetPool()));
    }

    // Compute (b + a·sk_small) mod q per prime, returned in coefficient form.
    // We use SEAL convention: decryption = c0 + c1·s.
    std::vector<std::vector<std::uint64_t>> r_per_prime(L, std::vector<std::uint64_t>(N, 0));
    for (std::size_t k = 0; k < L; k++) {
        const Modulus& q = rs.q_chain[k];
        // NTT-encode sk_small mod q_k.
        std::vector<std::uint64_t> sk_ntt(N);
        for (std::size_t i = 0; i < N; i++) {
            int v = sk_small_coeffs[i];
            sk_ntt[i] = (v >= 0) ? static_cast<std::uint64_t>(v)
                                 : (q.value() - static_cast<std::uint64_t>(-v));
        }
        ntt_negacyclic_harvey(sk_ntt.data(), *ntt_small[k]);

        // NTT-encode a (already coeff form from ring_switch).
        std::vector<std::uint64_t> a_ntt(rs.a[k]);
        ntt_negacyclic_harvey(a_ntt.data(), *ntt_small[k]);

        // a · sk_small in NTT domain.
        std::vector<std::uint64_t> as(N);
        dyadic_product_coeffmod(a_ntt.data(), sk_ntt.data(), N, q, as.data());

        // NTT-encode b.
        std::vector<std::uint64_t> b_ntt(rs.b[k]);
        ntt_negacyclic_harvey(b_ntt.data(), *ntt_small[k]);

        // r = b + a·sk_small (SEAL convention).
        std::vector<std::uint64_t> r_ntt(N);
        add_poly_coeffmod(b_ntt.data(), as.data(), N, q, r_ntt.data());

        // Inverse NTT to coefficient form mod q_k.
        inverse_ntt_negacyclic_harvey(r_ntt.data(), *ntt_small[k]);
        r_per_prime[k] = std::move(r_ntt);
    }

    // CRT-reconstruct each coefficient into a signed value mod Q, then reduce mod p.
    Modulus p_mod(plain_modulus);
    std::vector<std::uint64_t> mu(N, 0);

    if (L == 1) {
        const std::uint64_t q0 = rs.q_chain[0].value();
        const std::uint64_t q_half = q0 >> 1;
        const std::uint64_t q_mod_p = barrett_reduce_64(q0, p_mod);
        for (std::size_t i = 0; i < N; i++) {
            const std::uint64_t v = r_per_prime[0][i];
            const std::uint64_t v_mod_p = barrett_reduce_64(v, p_mod);
            if (v > q_half) {
                // negative branch: signed v - q ≡ v_mod_p - q_mod_p (mod p)
                mu[i] = (v_mod_p >= q_mod_p) ? (v_mod_p - q_mod_p)
                                              : (plain_modulus - (q_mod_p - v_mod_p));
            } else {
                mu[i] = v_mod_p;
            }
        }
    } else {
        // L == 2: uint128 CRT.
        const std::uint64_t q0 = rs.q_chain[0].value();
        const std::uint64_t q1 = rs.q_chain[1].value();
        const __uint128_t Q = static_cast<__uint128_t>(q0) * static_cast<__uint128_t>(q1);
        const __uint128_t Q_half = Q >> 1;
        std::uint64_t inv_q0_mod_q1 = 0;
        if (!try_invert_uint_mod(q0, rs.q_chain[1], inv_q0_mod_q1)) {
            throw std::runtime_error("decrypt_ring_switched: q0 not invertible mod q1");
        }
        const std::int64_t p_val = static_cast<std::int64_t>(plain_modulus);
        for (std::size_t i = 0; i < N; i++) {
            const std::uint64_t r0 = r_per_prime[0][i];
            const std::uint64_t r1 = r_per_prime[1][i];
            const std::uint64_t r0_mod_q1 = barrett_reduce_64(r0, rs.q_chain[1]);
            std::uint64_t y = (r1 >= r0_mod_q1) ? (r1 - r0_mod_q1) : (q1 - (r0_mod_q1 - r1));
            y = multiply_uint_mod(y, inv_q0_mod_q1, rs.q_chain[1]);
            __uint128_t x = static_cast<__uint128_t>(r0) + static_cast<__uint128_t>(q0) * static_cast<__uint128_t>(y);
            __int128_t xs = (x > Q_half) ? static_cast<__int128_t>(x) - static_cast<__int128_t>(Q)
                                          : static_cast<__int128_t>(x);
            std::int64_t r = static_cast<std::int64_t>(xs % p_val);
            if (r < 0) r += p_val;
            mu[i] = static_cast<std::uint64_t>(r);
        }
    }

    // The trace doubled the encoded plaintext: actual encoded value is 2·μ_orig · correction.
    // Mod-switches accumulated `correction = ∏ q_dropped (mod p)` in SEAL's convention; divide it out.
    // Combined inverse to apply: (2 · correction)⁻¹ mod p.
    std::uint64_t two_times_corr = multiply_uint_mod(2ULL, rs.correction, p_mod);
    std::uint64_t inv_factor = 0;
    if (!try_invert_uint_mod(two_times_corr, p_mod, inv_factor)) {
        throw std::runtime_error("decrypt_ring_switched: (2·correction) not invertible mod p");
    }
    multiply_poly_scalar_coeffmod(mu.data(), N, inv_factor, p_mod, mu.data());

    return mu;
}

std::vector<std::uint64_t> decrypt_and_measure_noise(const RingSwitchedCt& rs,
                                                     const std::vector<int>& sk_small_coeffs,
                                                     std::uint64_t plain_modulus,
                                                     double* out_log2_noise) {
    const std::size_t N = rs.N_small;
    const std::size_t L = rs.q_chain.size();
    if (sk_small_coeffs.size() != N) throw std::runtime_error("decrypt_and_measure: sk size mismatch");
    if (L != 1 && L != 2) throw std::runtime_error("decrypt_and_measure: only L in {1,2}");

    int log_N = 0; while ((std::size_t(1) << log_N) < N) log_N++;
    std::vector<std::unique_ptr<NTTTables>> ntt_small;
    ntt_small.reserve(L);
    for (std::size_t k = 0; k < L; k++) {
        ntt_small.emplace_back(std::make_unique<NTTTables>(log_N, rs.q_chain[k],
                                                            MemoryManager::GetPool()));
    }

    // Compute r = b + a·sk_small per prime, in coefficient form.
    std::vector<std::vector<std::uint64_t>> r_per_prime(L, std::vector<std::uint64_t>(N, 0));
    for (std::size_t k = 0; k < L; k++) {
        const Modulus& q = rs.q_chain[k];
        std::vector<std::uint64_t> sk_ntt(N);
        for (std::size_t i = 0; i < N; i++) {
            int v = sk_small_coeffs[i];
            sk_ntt[i] = (v >= 0) ? static_cast<std::uint64_t>(v)
                                 : (q.value() - static_cast<std::uint64_t>(-v));
        }
        ntt_negacyclic_harvey(sk_ntt.data(), *ntt_small[k]);
        std::vector<std::uint64_t> a_ntt(rs.a[k]);
        ntt_negacyclic_harvey(a_ntt.data(), *ntt_small[k]);
        std::vector<std::uint64_t> as(N);
        dyadic_product_coeffmod(a_ntt.data(), sk_ntt.data(), N, q, as.data());
        std::vector<std::uint64_t> b_ntt(rs.b[k]);
        ntt_negacyclic_harvey(b_ntt.data(), *ntt_small[k]);
        std::vector<std::uint64_t> r_ntt(N);
        add_poly_coeffmod(b_ntt.data(), as.data(), N, q, r_ntt.data());
        inverse_ntt_negacyclic_harvey(r_ntt.data(), *ntt_small[k]);
        r_per_prime[k] = std::move(r_ntt);
    }

    // CRT-reconstruct each coefficient into a SIGNED __int128 value (works for L ≤ 2 with q ≤ ~127 b).
    Modulus p_mod(plain_modulus);
    std::vector<std::uint64_t> mu_pre_unscale(N);    // 2·μ·correction mod p (before /2·corr)
    std::vector<__int128_t> v_signed(N);             // signed reconstruction mod Q

    if (L == 1) {
        const std::uint64_t q0 = rs.q_chain[0].value();
        const std::uint64_t q_half = q0 >> 1;
        for (std::size_t i = 0; i < N; i++) {
            const std::uint64_t v = r_per_prime[0][i];
            v_signed[i] = (v > q_half)
                ? static_cast<__int128_t>(v) - static_cast<__int128_t>(q0)
                : static_cast<__int128_t>(v);
            std::uint64_t v_mod_p = barrett_reduce_64(v, p_mod);
            std::uint64_t q_mod_p = barrett_reduce_64(q0, p_mod);
            if (v > q_half) {
                mu_pre_unscale[i] = (v_mod_p >= q_mod_p) ? (v_mod_p - q_mod_p)
                                                          : (plain_modulus - (q_mod_p - v_mod_p));
            } else {
                mu_pre_unscale[i] = v_mod_p;
            }
        }
    } else {
        const std::uint64_t q0 = rs.q_chain[0].value();
        const std::uint64_t q1 = rs.q_chain[1].value();
        const __uint128_t Q = static_cast<__uint128_t>(q0) * static_cast<__uint128_t>(q1);
        const __uint128_t Q_half = Q >> 1;
        std::uint64_t inv_q0_mod_q1 = 0;
        if (!try_invert_uint_mod(q0, rs.q_chain[1], inv_q0_mod_q1))
            throw std::runtime_error("decrypt_and_measure: q0 not invertible mod q1");
        const std::int64_t p_val_i = static_cast<std::int64_t>(plain_modulus);
        for (std::size_t i = 0; i < N; i++) {
            const std::uint64_t r0 = r_per_prime[0][i];
            const std::uint64_t r1 = r_per_prime[1][i];
            const std::uint64_t r0_mod_q1 = barrett_reduce_64(r0, rs.q_chain[1]);
            std::uint64_t y = (r1 >= r0_mod_q1) ? (r1 - r0_mod_q1) : (q1 - (r0_mod_q1 - r1));
            y = multiply_uint_mod(y, inv_q0_mod_q1, rs.q_chain[1]);
            __uint128_t x = static_cast<__uint128_t>(r0)
                          + static_cast<__uint128_t>(q0) * static_cast<__uint128_t>(y);
            __int128_t xs = (x > Q_half) ? static_cast<__int128_t>(x) - static_cast<__int128_t>(Q)
                                          : static_cast<__int128_t>(x);
            v_signed[i] = xs;
            std::int64_t r = static_cast<std::int64_t>(xs % p_val_i);
            if (r < 0) r += p_val_i;
            mu_pre_unscale[i] = static_cast<std::uint64_t>(r);
        }
    }

    // Recover μ = (mu_pre_unscale) · (2·correction)⁻¹ mod p.
    std::uint64_t two_times_corr = multiply_uint_mod(2ULL, rs.correction, p_mod);
    std::uint64_t inv_factor = 0;
    if (!try_invert_uint_mod(two_times_corr, p_mod, inv_factor))
        throw std::runtime_error("decrypt_and_measure: (2·correction) not invertible mod p");
    std::vector<std::uint64_t> mu(N);
    multiply_poly_scalar_coeffmod(mu_pre_unscale.data(), N, inv_factor, p_mod, mu.data());

    // Self-measure noise inf-norm using the just-decrypted μ as the reference.
    //   v_signed[i] = p·e[i] + (2·μ·correction mod p)_lifted   (mod Q, centered).
    //   Recover e[i] by trying both lifts of (2·μ_decoded[i]·correction mod p) into Z.
    if (out_log2_noise) {
        __int128_t max_abs = 0;
        const std::int64_t p_val_i = static_cast<std::int64_t>(plain_modulus);
        for (std::size_t i = 0; i < N; i++) {
            std::uint64_t enc = multiply_uint_mod(mu[i], two_times_corr, p_mod);
            __int128_t diff = v_signed[i] - static_cast<__int128_t>(enc);
            if (diff % p_val_i != 0) {
                __int128_t alt = v_signed[i] - (static_cast<__int128_t>(enc) - p_val_i);
                if (alt % p_val_i == 0) diff = alt;
                else { *out_log2_noise = -1.0; return mu; }
            }
            __int128_t e = diff / p_val_i;
            if (e < 0) e = -e;
            if (e > max_abs) max_abs = e;
        }
        if (max_abs == 0) {
            *out_log2_noise = 0.0;
        } else {
            int hi = 0;
            __uint128_t u = static_cast<__uint128_t>(max_abs);
            while (u > 1) { u >>= 1; hi++; }
            double frac = std::log2(static_cast<double>(static_cast<__uint128_t>(max_abs)
                                                       / (static_cast<__uint128_t>(1) << hi)));
            *out_log2_noise = hi + frac;
        }
    }

    return mu;
}

}  // namespace ring_switch_util
