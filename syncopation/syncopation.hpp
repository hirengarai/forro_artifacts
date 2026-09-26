/*
 * syncopation.hpp
 *
 * Algorithm 1 of Wang, Liu, Hou, Lin, CRYPTO 2023
 * (https://eprint.iacr.org/2023/1087), for Forro: given a PNB set, the segments
 * of the backward computation that become writable once a restriction holds.
 *
 * by Hiren
 *
 * Not a program: included by verify.hpp. Build with  make  in this folder.
 */

#pragma once

#include <array>
#include <bit>
#include <cstddef>
#include <vector>

// u8/u16/u32 and the Forro layout come from forro.cpp, which verify.hpp
// includes before this file.

namespace synco
{
    inline constexpr size_t STATE_WORDS = 16;
    inline constexpr size_t WORD_BITS   = 32;

    // Forro's key sits in x[0..3] and x[8..11]. Only T0 applies: the deeper
    // levels need an unrotated subtraction, which Forro's inverse SRF lacks.
    inline constexpr std::array<u8, 8> FORRO_LAYOUT = {0, 1, 2, 3, 8, 9, 10, 11};

    struct Condition
    {
        u8  word;
        u8  restricted_bit;
        u8  index;
        u32 seg_mask;
        u32 carry_mask;
    };

    struct Plan
    {
        std::array<u32, STATE_WORDS> seg{};
        std::array<u32, STATE_WORDS> carry{};
        std::array<u32, STATE_WORDS> free{};   // writable with no condition
        std::vector<Condition>       conditions;
        u32                          restricted_bits = 0;

        size_t theta() const { return conditions.size(); }
    };

    inline constexpr u32 range_mask(u32 lo, u32 hi)
    {
        if (hi < lo)
            return 0;
        const u32 width = hi - lo + 1;
        return (width >= WORD_BITS ? 0xFFFFFFFFu : ((1u << width) - 1u)) << lo;
    }

    // Per key word, walk the runs of non-PNB bits. A run above a PNB spends its
    // lowest bit as the restriction and yields the segment above it plus the
    // carry bit on top. A run at bit 0 is free. A run whose segment is empty is
    // still a condition when it buys the carry bit.
    inline Plan build_plan(const std::array<u32, STATE_WORDS> &pnb)
    {
        Plan plan;
        u8   index = 0;

        for (const u8 w : FORRO_LAYOUT)
        {
            const u32 pnb_mask = pnb[w];

            u32 j = 0;
            while (j < WORD_BITS)
            {
                if ((pnb_mask >> j) & 1u)
                {
                    ++j;
                    continue;
                }

                const u32 run_lo = j;
                while (j < WORD_BITS && !((pnb_mask >> j) & 1u))
                    ++j;
                const u32 run_hi = j - 1;

                if (run_lo == 0)
                {
                    plan.seg[w]  |= range_mask(0, run_hi);
                    plan.free[w] |= range_mask(0, run_hi);
                    continue;
                }

                const u32 seg   = range_mask(run_lo + 1, run_hi);
                const u32 carry = (run_hi + 1 < WORD_BITS)
                                      ? range_mask(run_hi + 1, run_hi + 1) : 0u;
                if (!seg && !carry)
                    continue;

                plan.restricted_bits += 1;
                plan.seg[w]   |= seg;
                plan.carry[w] |= carry;
                plan.conditions.push_back({w, static_cast<u8>(run_lo), index, seg, carry});
                ++index;
            }
        }
        return plan;
    }

    // Rebuild seg/carry from the unconditional part plus the conditions still
    // active, and recount |k^R|. Used by pruning.
    inline void apply_active(Plan &plan, const std::vector<char> &active)
    {
        plan.seg   = plan.free;
        plan.carry = {};

        std::array<u32, STATE_WORDS> restricted{};
        for (size_t c = 0; c < plan.conditions.size(); ++c)
        {
            if (!active[c])
                continue;
            const Condition &cond = plan.conditions[c];
            plan.seg[cond.word]   |= cond.seg_mask;
            plan.carry[cond.word] |= cond.carry_mask;
            restricted[cond.word] |= (1u << cond.restricted_bit);
        }

        plan.restricted_bits = 0;
        for (u32 r : restricted)
            plan.restricted_bits += static_cast<u32>(std::popcount(r));
    }

    struct Assignment
    {
        std::array<u32, STATE_WORDS> zero{};
        std::array<u32, STATE_WORDS> one{};
    };

    // g must be a fixed function for Proposition 1 to hold.
    //
    // pattern: each block of consecutive PNBs gets 100...0 (top bit 1, rest 0),
    //          singletons get 0. This minimises the probability that the bit
    //          above the block is disturbed -- 1/4 for every block size, against
    //          1/2 - 2^-(m+1) for all-zero. Optimal: exhausting all 2^m patterns
    //          gives 1/4 only for 100...0 and its complement 011...1.
    // zero   : every PNB gets 0, the convention of Aumasson et al. Measured 2^3
    //          to 2^4 worse in eps_a. Only for reproducing published numbers.
    template <typename ToWordBit>
    inline Assignment assign(const std::vector<u16> &sorted_set, ToWordBit to_word_bit,
                             bool use_pattern_split)
    {
        Assignment a;

        auto set = [&](u16 idx, std::array<u32, STATE_WORDS> &dst)
        {
            u16 w = 0, b = 0;
            to_word_bit(idx, w, b);
            dst[w] |= (1u << b);
        };

        if (!use_pattern_split)
        {
            for (u16 idx : sorted_set)
                set(idx, a.zero);
            return a;
        }

        for (size_t i = 0; i < sorted_set.size();)
        {
            // A run must stay inside one 32-bit key word: the pattern is about a
            // carry chain, and carries do not cross from k_w[31] into k_{w+1}[0].
            // Key bits 191 and 192 are adjacent as indices but sit in different
            // words, so testing index adjacency alone merges unrelated runs.
            size_t j = i;
            while (j + 1 < sorted_set.size()
                   && sorted_set[j + 1] == sorted_set[j] + 1
                   && sorted_set[j + 1] / WORD_BITS == sorted_set[j] / WORD_BITS)
                ++j;

            if (j > i)
            {
                for (size_t t = i; t < j; ++t)
                    set(sorted_set[t], a.zero);
                set(sorted_set[j], a.one);
            }
            else
                set(sorted_set[i], a.zero);

            i = j + 1;
        }
        return a;
    }
}
