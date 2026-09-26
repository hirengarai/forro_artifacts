/*
 * forro.cpp
 *
 * Forro, self-contained: the fixed-width aliases, a small u32 SIMD layer, the
 * cipher forward and backward, and a self-check. Nothing outside this file is
 * needed by the two verification programs, which include it directly.
 *
 * Forro (Coutinho et al., ASIACRYPT 2022) has a 4x4 state of 32-bit words:
 *
 *      k0  k1  k2  k3          key      : x[0..3], x[8..11]
 *      t0  t1  c0  c1          counter  : x[4], x[5], x[12], x[13]
 *      k8  k9 k10 k11          nonce    : the other two of those four
 *      t2  t3  c2  c3          constants: x[6], x[7], x[14], x[15]
 *
 * The round function applies the subround function (SRF) to four tuples of five
 * words. The fifth word is the Pollination word. One tuple is therefore 0.25
 * round, and that is the atomic step here: the attacks invert 10 or 11
 * subrounds, not whole rounds, so every index below counts tuples.
 *
 * Odd rounds use the column tuples, even rounds the diagonal tuples. The SRF is
 * three ARX steps, rotating by 10, 27 and 8:
 *
 *      ARX10:  d += e;  c ^= d;  b = rotl(b + c, 10)
 *      ARX27:  a += b;  e ^= a;  d = rotl(d + e, 27)
 *      ARX8 :  c += d;  b ^= c;  a = rotl(a + b,  8)
 *
 * Monte-Carlo sampling uses one SIMD lane per independent sample -- the sixteen
 * state words become sixteen vectors and every step is lane-wise, so lanes never
 * interact. The backend is chosen automatically: NEON (4 lanes), AVX2 (8 lanes),
 * or a scalar fallback (1 lane, or forced with -DFORCE_SCALAR).
 *
 * by Hiren
 *
 * Not a program: included by forro65_verification.cpp and forro675_verification.cpp.
 */

#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <initializer_list>

using u8  = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

// ============================== SIMD layer ==============================
//
// Only the operations the attacks use. Rotations are templated on the amount
// because every ISA needs a compile-time literal, which is why the cipher steps
// below are templated on the tuple index rather than dispatching at runtime.

// ------- NEON -------

#if !defined(FORCE_SCALAR) && (defined(__ARM_NEON) || defined(__ARM_NEON__))
#include <arm_neon.h>
#define SIMD_ISA "NEON"

namespace simd
{
    using vec = uint32x4_t;
    constexpr std::size_t LANES = 4;

    inline vec set1(u32 v)        { return vdupq_n_u32(v); }
    inline vec add(vec a, vec b)  { return vaddq_u32(a, b); }
    inline vec sub(vec a, vec b)  { return vsubq_u32(a, b); }
    inline vec bxor(vec a, vec b) { return veorq_u32(a, b); }
    inline vec band(vec a, vec b) { return vandq_u32(a, b); }
    inline vec bor(vec a, vec b)  { return vorrq_u32(a, b); }
    inline vec mul(vec a, vec b)  { return vmulq_u32(a, b); }
    inline vec andnot(vec a, vec b) { return vbicq_u32(a, b); }   // a & ~b

    template <int N> inline vec shl(vec a) { return vshlq_n_u32(a, N); }

    // VSRI keeps the upper (32-N) bits of its first operand -- exactly the
    // (x << N) half -- and inserts (x >> (32-N)) underneath. Two ops.
    template <int N> inline vec rotl(vec x)
    {
        static_assert(N > 0 && N < 32, "rotate amount must be in [1,31]");
        return vsriq_n_u32(vshlq_n_u32(x, N), x, 32 - N);
    }

    // Per-lane parity of the full 32-bit word, returned in bit 0.
    inline vec parity(vec t)
    {
        const uint8x16_t bytes = vcntq_u8(vreinterpretq_u8_u32(t));
        return vandq_u32(vpaddlq_u16(vpaddlq_u8(bytes)), vdupq_n_u32(1));
    }

    inline u32  hsum(vec a)            { return vaddvq_u32(a); }
    inline void store(u32 *dst, vec a) { vst1q_u32(dst, a); }
    inline vec  load(const u32 *src)   { return vld1q_u32(src); }
}

// ------- AVX2 -------

#elif !defined(FORCE_SCALAR) && defined(__AVX2__)
#include <immintrin.h>
#define SIMD_ISA "AVX2"

namespace simd
{
    using vec = __m256i;
    constexpr std::size_t LANES = 8;

    inline vec set1(u32 v)        { return _mm256_set1_epi32(static_cast<int>(v)); }
    inline vec add(vec a, vec b)  { return _mm256_add_epi32(a, b); }
    inline vec sub(vec a, vec b)  { return _mm256_sub_epi32(a, b); }
    inline vec bxor(vec a, vec b) { return _mm256_xor_si256(a, b); }
    inline vec band(vec a, vec b) { return _mm256_and_si256(a, b); }
    inline vec bor(vec a, vec b)  { return _mm256_or_si256(a, b); }
    inline vec mul(vec a, vec b)  { return _mm256_mullo_epi32(a, b); }
    inline vec andnot(vec a, vec b) { return _mm256_andnot_si256(b, a); }   // a & ~b

    template <int N> inline vec shl(vec a) { return _mm256_slli_epi32(a, N); }

    template <int N> inline vec rotl(vec x)
    {
        static_assert(N > 0 && N < 32, "rotate amount must be in [1,31]");
        // 16 and 8 are pure byte permutations -- one op instead of three.
        if constexpr (N == 16)
            return _mm256_shuffle_epi8(
                x, _mm256_setr_epi8(2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13,
                                    2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13));
        else if constexpr (N == 8)
            return _mm256_shuffle_epi8(
                x, _mm256_setr_epi8(3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10, 15, 12, 13, 14,
                                    3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10, 15, 12, 13, 14));
        else
            return _mm256_or_si256(_mm256_slli_epi32(x, N), _mm256_srli_epi32(x, 32 - N));
    }

    // Per-lane parity in bit 0, via the XOR-fold (AVX2 has no per-lane popcount).
    inline vec parity(vec t)
    {
        t = _mm256_xor_si256(t, _mm256_srli_epi32(t, 16));
        t = _mm256_xor_si256(t, _mm256_srli_epi32(t, 8));
        t = _mm256_xor_si256(t, _mm256_srli_epi32(t, 4));
        t = _mm256_xor_si256(t, _mm256_srli_epi32(t, 2));
        t = _mm256_xor_si256(t, _mm256_srli_epi32(t, 1));
        return _mm256_and_si256(t, _mm256_set1_epi32(1));
    }

    inline u32 hsum(vec a)
    {
        alignas(32) u32 tmp[LANES];
        _mm256_store_si256(reinterpret_cast<__m256i *>(tmp), a);
        u32 s = 0;
        for (std::size_t i = 0; i < LANES; ++i)
            s += tmp[i];
        return s;
    }

    inline void store(u32 *dst, vec a)
    {
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(dst), a);
    }
    inline vec load(const u32 *src)
    {
        return _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src));
    }
}

// ------- scalar fallback -------

#else
#define SIMD_ISA "scalar"
// simd::vec is u32 here, so any function overloaded on both scalars and vectors
// would collide. Guard those with this.
#define SIMD_VEC_IS_U32 1

// Landing here without asking for it is almost always a missing compiler flag,
// and it costs the whole SIMD speedup silently -- the build still succeeds and
// the answers are still correct, just 4-8x slower. Say so loudly. Pass
// -DFORCE_SCALAR when the 1-lane build is what you actually want.
#if !defined(FORCE_SCALAR)
#warning "forro.cpp: no SIMD backend detected, building the 1-lane scalar fallback. \
On x86-64 add -mavx2 (or -march=native); on aarch64 check -mcpu=native. \
Pass -DFORCE_SCALAR to silence this."
#endif

namespace simd
{
    using vec = u32;
    constexpr std::size_t LANES = 1;

    inline vec set1(u32 v)        { return v; }
    inline vec add(vec a, vec b)  { return static_cast<u32>(a + b); }
    inline vec sub(vec a, vec b)  { return static_cast<u32>(a - b); }
    inline vec bxor(vec a, vec b) { return a ^ b; }
    inline vec band(vec a, vec b) { return a & b; }
    inline vec bor(vec a, vec b)  { return a | b; }
    inline vec mul(vec a, vec b)  { return static_cast<u32>(a * b); }
    inline vec andnot(vec a, vec b) { return a & ~b; }

    template <int N> inline vec shl(vec a) { return static_cast<u32>(a << N); }

    template <int N> inline vec rotl(vec x)
    {
        static_assert(N > 0 && N < 32, "rotate amount must be in [1,31]");
        return std::rotl(x, N);
    }

    inline vec  parity(vec t)          { return static_cast<u32>(std::popcount(t) & 1); }
    inline u32  hsum(vec a)            { return a; }
    inline void store(u32 *dst, vec a) { *dst = a; }
    inline vec  load(const u32 *src)   { return *src; }
}
#endif

namespace simd
{
    template <int N> inline vec rotr(vec x) { return rotl<32 - N>(x); }

    // xoshiro128** (https://prng.di.unimi.it/xoshiro128starstar.c) with one
    // independent stream per lane, seeded through SplitMix64. A per-call
    // std::uniform_int_distribution is a no-op transform over a full-range draw,
    // but it is not always inlined away and becomes a large share of the cost
    // once the arithmetic core is vectorised.
    struct Xoshiro128StarStar
    {
        vec s0, s1, s2, s3;

        explicit Xoshiro128StarStar(u64 seed)
        {
            alignas(64) u32 w[4][LANES];
            u64 sm = seed;

            auto splitmix = [&sm]() -> u64
            {
                sm += 0x9e3779b97f4a7c15ULL;
                u64 z = sm;
                z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
                z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
                return z ^ (z >> 31);
            };

            for (std::size_t l = 0; l < LANES; ++l)
                for (int i = 0; i < 4; ++i)
                {
                    const u32 v = static_cast<u32>(splitmix() >> 32);
                    w[i][l] = v ? v : 0x9e3779b9u;   // state must not be all-zero
                }

            s0 = load(w[0]);
            s1 = load(w[1]);
            s2 = load(w[2]);
            s3 = load(w[3]);
        }

        inline vec next()
        {
            const vec result = mul(rotl<7>(mul(s1, set1(5))), set1(9));
            const vec t      = shl<9>(s1);

            s2 = bxor(s2, s0);
            s3 = bxor(s3, s1);
            s1 = bxor(s1, s2);
            s0 = bxor(s0, s3);
            s2 = bxor(s2, t);
            s3 = rotl<11>(s3);

            return result;
        }
    };
}

// ============================== the cipher ==============================

namespace forro
{
    constexpr std::size_t STATE_WORDS      = 16;
    constexpr std::size_t KEY_COUNT        = 8;
    constexpr std::size_t KEY_LO_START     = 0;    // x[0..3]
    constexpr std::size_t KEY_HI_START     = 8;    // x[8..11]
    constexpr int         TUPLES_PER_ROUND = 4;    // one tuple = 0.25 round

    constexpr const char *NAME = "Forro";

    // Constants at x[6], x[7], x[14], x[15]: "voltamos a sabranca" in ASCII.
    constexpr u32 C0 = 0x746C6F76;
    constexpr u32 C1 = 0x61616461;
    constexpr u32 C2 = 0x72626173;
    constexpr u32 C3 = 0x61636E61;

    constexpr int ROT_A = 10;   // ARX10
    constexpr int ROT_B = 27;   // ARX27
    constexpr int ROT_C = 8;    // ARX8

    // The five word indices (a, b, c, d, e) of each SRF tuple. e is the
    // Pollination word, and it is shared with a neighbouring tuple -- which is
    // exactly what makes the trail bookkeeping non-routine.
    inline constexpr int COL[4][5]  = {{0, 4,  8, 12, 3}, {1, 5,  9, 13, 0},
                                       {2, 6, 10, 14, 1}, {3, 7, 11, 15, 2}};
    inline constexpr int DIAG[4][5] = {{0, 5, 10, 15, 3}, {1, 6, 11, 12, 0},
                                       {2, 7,  8, 13, 1}, {3, 4,  9, 14, 2}};

    // Feed-forward: the keystream is Z = X^(R) + X^(0), word-wise.
    inline u32 ff(u32 fin, u32 ini)     { return static_cast<u32>(fin + ini); }
    inline u32 ff_inv(u32 z, u32 ini)   { return static_cast<u32>(z - ini); }

#ifndef SIMD_VEC_IS_U32
    inline simd::vec ff(simd::vec fin, simd::vec ini)   { return simd::add(fin, ini); }
    inline simd::vec ff_inv(simd::vec z, simd::vec ini) { return simd::sub(z, ini); }
#endif

    // Key bit index -> (state word, bit). Forro keeps its key in the two halves
    // x[0..3] and x[8..11], so bit 128 lands in x[8], not x[4].
    inline void key_bit_to_word(u16 index, u16 &word, u16 &bit)
    {
        const u16 key_word = index / 32;
        bit = index % 32;
        word = (index < 128) ? u16(key_word + KEY_LO_START)
                             : u16(key_word - 4 + KEY_HI_START);
    }

    // Inverse of key_bit_to_word.
    inline u16 word_bit_to_key_bit(u8 word, u8 bit)
    {
        const u16 kw = (word < KEY_HI_START) ? u16(word) : u16(word - KEY_HI_START + 4);
        return u16(kw * 32 + bit);
    }

    inline void copy_state(simd::vec *dst, const simd::vec *src)
    {
        std::memcpy(dst, src, STATE_WORDS * sizeof(simd::vec));
    }

    namespace detail
    {
        // The five word indices of tuple T in the odd (COL) or even (DIAG)
        // layer, resolved entirely at compile time.
        template <bool ODD, int T> struct Tuple
        {
            static constexpr int a = ODD ? COL[T][0] : DIAG[T][0];
            static constexpr int b = ODD ? COL[T][1] : DIAG[T][1];
            static constexpr int c = ODD ? COL[T][2] : DIAG[T][2];
            static constexpr int d = ODD ? COL[T][3] : DIAG[T][3];
            static constexpr int e = ODD ? COL[T][4] : DIAG[T][4];
        };

        // true when tuple index K falls in an odd (COL) round.
        template <int K>
        inline constexpr bool is_odd_round = ((K / TUPLES_PER_ROUND + 1) & 1) != 0;
    }

    // ------- SIMD forward -------

    struct VecForward
    {
        template <bool ODD, int T>
        static inline void srf(simd::vec *x)
        {
            using L = detail::Tuple<ODD, T>;

            // ARX10
            x[L::d] = simd::add(x[L::d], x[L::e]);
            x[L::c] = simd::bxor(x[L::c], x[L::d]);
            x[L::b] = simd::rotl<ROT_A>(simd::add(x[L::b], x[L::c]));

            // ARX27
            x[L::a] = simd::add(x[L::a], x[L::b]);
            x[L::e] = simd::bxor(x[L::e], x[L::a]);
            x[L::d] = simd::rotl<ROT_B>(simd::add(x[L::d], x[L::e]));

            // ARX8
            x[L::c] = simd::add(x[L::c], x[L::d]);
            x[L::b] = simd::bxor(x[L::b], x[L::c]);
            x[L::a] = simd::rotl<ROT_C>(simd::add(x[L::a], x[L::b]));
        }

        // One atomic 0.25-round tick: the SRF on tuple (K % 4) of the layer that
        // round K/4+1 selects.
        template <int K>
        static inline void step(simd::vec *x)
        {
            srf<detail::is_odd_round<K>, K % TUPLES_PER_ROUND>(x);
        }

        // Tuples [FROM, TO), fully unrolled.
        template <int FROM, int TO>
        static inline void range(simd::vec *x)
        {
            if constexpr (FROM < TO)
            {
                step<FROM>(x);
                range<FROM + 1, TO>(x);
            }
        }
    };

    // ------- SIMD backward -------

    struct VecBackward
    {
        // Exact inverse of VecForward::srf<ODD, T>, undoing ARX8, ARX27, ARX10
        // in that order.
        template <bool ODD, int T>
        static inline void srf(simd::vec *x)
        {
            using L = detail::Tuple<ODD, T>;

            // inverse ARX8
            x[L::a] = simd::sub(simd::rotr<ROT_C>(x[L::a]), x[L::b]);
            x[L::b] = simd::bxor(x[L::b], x[L::c]);
            x[L::c] = simd::sub(x[L::c], x[L::d]);

            // inverse ARX27
            x[L::d] = simd::sub(simd::rotr<ROT_B>(x[L::d]), x[L::e]);
            x[L::e] = simd::bxor(x[L::e], x[L::a]);
            x[L::a] = simd::sub(x[L::a], x[L::b]);

            // inverse ARX10
            x[L::b] = simd::sub(simd::rotr<ROT_A>(x[L::b]), x[L::c]);
            x[L::c] = simd::bxor(x[L::c], x[L::d]);
            x[L::d] = simd::sub(x[L::d], x[L::e]);
        }

        template <int K>
        static inline void step(simd::vec *x)
        {
            srf<detail::is_odd_round<K>, K % TUPLES_PER_ROUND>(x);
        }

        // Undoes tuples [FROM, TO) in reverse order, fully unrolled.
        template <int FROM, int TO>
        static inline void range(simd::vec *x)
        {
            if constexpr (FROM < TO)
            {
                step<TO - 1>(x);
                range<FROM, TO - 1>(x);
            }
        }
    };

    // ------- scalar reference -------
    //
    // Used only by self_check(), as an independent second implementation to hold
    // the SIMD path against. It indexes COL/DIAG at run time, so it shares no
    // code with the templated version above.

    inline void srf_scalar(u32 &a, u32 &b, u32 &c, u32 &d, u32 &e)
    {
        d += e; c ^= d; b = std::rotl(u32(b + c), ROT_A);
        a += b; e ^= a; d = std::rotl(u32(d + e), ROT_B);
        c += d; b ^= c; a = std::rotl(u32(a + b), ROT_C);
    }

    inline void srf_inv_scalar(u32 &a, u32 &b, u32 &c, u32 &d, u32 &e)
    {
        a = u32(std::rotr(a, ROT_C) - b); b ^= c; c -= d;
        d = u32(std::rotr(d, ROT_B) - e); e ^= a; a -= b;
        b = u32(std::rotr(b, ROT_A) - c); c ^= d; d -= e;
    }

    inline void step_scalar(u32 *x, int k, bool inverse)
    {
        const int t = k % TUPLES_PER_ROUND;
        const auto &idx = (((k / TUPLES_PER_ROUND + 1) & 1) != 0) ? COL[t] : DIAG[t];
        if (inverse)
            srf_inv_scalar(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]]);
        else
            srf_scalar(x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], x[idx[4]]);
    }

    inline void forward_scalar(u32 *x, int from, int to)
    {
        for (int k = from; k < to; ++k)
            step_scalar(x, k, false);
    }

    inline void backward_scalar(u32 *x, int from, int to)
    {
        for (int k = to - 1; k >= from; --k)
            step_scalar(x, k, true);
    }

    // ------- self-check -------

    // Two properties, on the test-vector state of the Forro specification:
    //
    //   1. backward undoes forward, over the full 14 rounds and over the
    //      fractional ranges the attacks actually invert;
    //   2. the SIMD path agrees with the scalar reference, lane by lane.
    //
    // Everything in the attacks rests on these two, so both verification
    // programs run this before measuring anything.
    inline bool self_check(bool verbose = false)
    {
        constexpr int FULL = 14 * TUPLES_PER_ROUND;   // 14 rounds = 56 tuples

        // "minha vida e dan / do par a o forro / e ste paisrbasanca"
        const u32 seed_state[STATE_WORDS] = {
            0x686E696D, 0x69762061, 0x65206164, 0x646E6120,
            0x00000000, 0x00000000, C0,         C1,
            0x70207261, 0x6520726F, 0x20657473, 0x73696170,
            0x74736F6D, 0x61206F72, C2,         C3};

        // 1. scalar round trip over the full 14 rounds
        u32 s[STATE_WORDS];
        std::memcpy(s, seed_state, sizeof(s));
        forward_scalar(s, 0, FULL);
        u32 forward_result[STATE_WORDS];
        std::memcpy(forward_result, s, sizeof(s));
        backward_scalar(s, 0, FULL);
        for (std::size_t w = 0; w < STATE_WORDS; ++w)
            if (s[w] != seed_state[w])
            {
                std::printf("[FAIL] scalar round trip differs at x[%zu]\n", w);
                return false;
            }

        // 2. SIMD forward matches the scalar forward in every lane
        simd::vec v[STATE_WORDS];
        for (std::size_t w = 0; w < STATE_WORDS; ++w)
            v[w] = simd::set1(seed_state[w]);
        VecForward::range<0, FULL>(v);
        for (std::size_t w = 0; w < STATE_WORDS; ++w)
        {
            u32 lanes[simd::LANES];
            simd::store(lanes, v[w]);
            for (std::size_t l = 0; l < simd::LANES; ++l)
                if (lanes[l] != forward_result[w])
                {
                    std::printf("[FAIL] SIMD lane %zu differs from scalar at x[%zu]\n", l, w);
                    return false;
                }
        }

        // 3. SIMD backward undoes it, over the two fractional ranges the attacks
        //    invert: 26 tuples (6.5 rounds) and 27 (6.75), back to tuple 16.
        for (const int total : {26, 27})
        {
            for (std::size_t w = 0; w < STATE_WORDS; ++w)
                v[w] = simd::set1(seed_state[w]);
            if (total == 26)
            {
                VecForward::range<0, 26>(v);
                VecBackward::range<16, 26>(v);
            }
            else
            {
                VecForward::range<0, 27>(v);
                VecBackward::range<16, 27>(v);
            }
            // The state is now at tuple 16; run it forward again and compare.
            u32 got[STATE_WORDS];
            u32 want[STATE_WORDS];
            std::memcpy(want, seed_state, sizeof(want));
            forward_scalar(want, 0, 16);
            for (std::size_t w = 0; w < STATE_WORDS; ++w)
            {
                u32 lanes[simd::LANES];
                simd::store(lanes, v[w]);
                got[w] = lanes[0];
                if (got[w] != want[w])
                {
                    std::printf("[FAIL] %d-tuple round trip differs at x[%zu]\n", total, w);
                    return false;
                }
            }
        }

        if (verbose)
            std::printf("self-check: %s, %d rounds forward and backward, "
                        "SIMD == scalar on %zu lanes -- OK\n",
                        NAME, FULL / TUPLES_PER_ROUND, simd::LANES);
        return true;
    }
}
