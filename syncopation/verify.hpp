/*
 * verify.hpp
 *
 * Verification of the conditional backward correlation eps'_a of a CPNB set of
 * Forro. This is the shared core of the two programs;
 * forro65_verification.cpp and forro675_verification.cpp only set the
 * parameters below and call verify::run.
 *
 * The estimator follows Proposition 1 of Wang, Liu, Hou and Lin, "Moving a Step
 * of ChaCha in Syncopated Rhythm", CRYPTO 2023 (https://eprint.iacr.org/2023/1087):
 * eps_a = (eps'_a)^2, so eps'_a is measured on a SINGLE initialisation and no
 * input difference is drawn. Reported as the median over keys of |correlation|.
 *
 * by Hiren
 *
 * Not a program: included by the two *_verification.cpp files.
 * Build both with  make  in this folder.
 */

#pragma once

// ---------------- parameters set by the including .cpp ----------------

#ifndef FORRO_TOTAL_ROUNDS_X4
#error "define FORRO_TOTAL_ROUNDS_X4 (rounds x 4) before including verify.hpp"
#endif
#ifndef FORRO_CPNB_FILE
#error "define FORRO_CPNB_FILE before including verify.hpp"
#endif
#ifndef FORRO_DEFAULT_LOG2_KEYS
#define FORRO_DEFAULT_LOG2_KEYS 12
#endif
#ifndef FORRO_DEFAULT_LOG2_IVS
#define FORRO_DEFAULT_LOG2_IVS 18
#endif

#include "forro.cpp"
#include "syncopation.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <future>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace verify
{
    using vec    = simd::vec;
    using VecRng = simd::Xoshiro128StarStar;
    using BitPos = std::pair<u16, u16>;

    inline constexpr size_t LANES = simd::LANES;

    // ---------------- attack parameters ----------------

    inline constexpr size_t KEY_SIZE_BITS = 256;

    // The 4-round differential-linear distinguisher:
    //   Gamma_out = v2[0] + v9[0] + v14[0] + v14[27],   eps_d = 2^{-16.22}.
    // Both attacks use it; they differ only in how many subrounds are inverted.
    inline constexpr double DISTINGUISHING_ROUND = 4.0;
    inline constexpr double FORWARD_CORR_LOG2    = -16.22;

    inline const std::vector<BitPos> OUTPUT_MASK_BITS = {{2, 0}, {9, 0}, {14, 0}, {14, 27}};

    // How g pins the CPNBs. Each block of consecutive CPNBs gets 100...0; see
    // synco::assign for why that beats the all-zero convention.
    inline constexpr bool USE_PATTERN_SPLIT = true;

    // ----------------------------------------------------

    inline constexpr size_t STATE_WORDS      = forro::STATE_WORDS;
    inline constexpr size_t KEY_COUNT        = forro::KEY_COUNT;
    inline constexpr size_t NUM_KEY_BITS     = 256;
    inline constexpr int    TUPLES_PER_ROUND = forro::TUPLES_PER_ROUND;

    inline constexpr int    TOTAL_STEP   = FORRO_TOTAL_ROUNDS_X4;
    inline constexpr double TOTAL_ROUNDS = double(TOTAL_STEP) / TUPLES_PER_ROUND;
    inline constexpr int    DIST_STEP    = int(DISTINGUISHING_ROUND * TUPLES_PER_ROUND);

    static_assert(TOTAL_STEP > DIST_STEP, "rounds must exceed the distinguishing round");

    inline constexpr std::array<size_t, 4> IV_WORDS = {4, 5, 12, 13};

    inline std::array<u32, STATE_WORDS> g_mask{};

    inline void init_masks()
    {
        g_mask.fill(0);
        for (const auto &b : OUTPUT_MASK_BITS)
            g_mask[b.first] |= (1u << b.second);
    }

    struct Candidate
    {
        synco::Assignment            pin;
        std::array<u32, STATE_WORDS> flip{};
        synco::Plan                  plan;
    };

    // ---------------- kernel ----------------

    inline void init_state(vec *x, const u32 *key, VecRng &rng)
    {
        x[6]  = simd::set1(forro::C0);
        x[7]  = simd::set1(forro::C1);
        x[14] = simd::set1(forro::C2);
        x[15] = simd::set1(forro::C3);
        for (size_t i = 0; i < KEY_COUNT; ++i)
            x[synco::FORRO_LAYOUT[i]] = simd::set1(key[i]);
        for (const size_t w : IV_WORDS)
            x[w] = rng.next();
    }

    // Keys are derived from (seed, index) so a run is reproducible and does not
    // depend on the thread count.
    inline void derive_key(u64 seed, u64 index, u32 *key)
    {
        u64 state = seed ^ (0x9E3779B97F4A7C15ULL * (index + 1));
        for (size_t i = 0; i < KEY_COUNT; ++i)
        {
            state += 0x9E3779B97F4A7C15ULL;
            u64 z = state;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            key[i] = static_cast<u32>((z ^ (z >> 31)) >> 32);
        }
    }

    // Segment bits take the true value (Lemma 2). The carry bit above takes
    // truth ^ key_diff: f and g subtract from the same known s, so
    // s[t] = x[t] ^ k[t] ^ Carry[t], the restriction fixes Carry[t], and bit t is
    // itself a CPNB. Writing the truth there loses correlation silently.
    inline void write_segments(vec *g_state, const vec *truth, const vec *key_diff,
                               const synco::Plan &plan)
    {
        for (const u8 w : synco::FORRO_LAYOUT)
        {
            const u32 mask = plan.seg[w] | plan.carry[w];
            if (!mask)
                continue;
            const vec want =
                plan.carry[w]
                    ? simd::bxor(truth[w], simd::band(key_diff[w], simd::set1(plan.carry[w])))
                    : truth[w];
            g_state[w] = simd::bxor(g_state[w],
                                    simd::band(simd::bxor(g_state[w], want), simd::set1(mask)));
        }
    }

    inline vec masked_value(const vec *s)
    {
        vec fold = simd::set1(0);
        for (size_t w = 0; w < STATE_WORDS; ++w)
            if (g_mask[w])
                fold = simd::bxor(fold, simd::band(s[w], simd::set1(g_mask[w])));
        return fold;
    }

    // One initialisation. f is the backward computation with the true key, g the
    // same with the CPNBs pinned and the syncopated segments written.
    inline void run_branch(const vec *x0, const Candidate *cands, size_t n_cand,
                           vec &fold_f, vec *fold_g)
    {
        vec x_final[STATE_WORDS], keystream[STATE_WORDS];
        vec f_state[STATE_WORDS], g_state[STATE_WORDS];
        vec key_diff[STATE_WORDS];

        forro::copy_state(x_final, x0);
        forro::VecForward::range<0, TOTAL_STEP>(x_final);
        for (size_t w = 0; w < STATE_WORDS; ++w)
            keystream[w] = forro::ff(x_final[w], x0[w]);

        forro::copy_state(f_state, x_final);
        forro::VecBackward::range<DIST_STEP, TOTAL_STEP>(f_state);
        fold_f = masked_value(f_state);

        for (size_t i = 0; i < n_cand; ++i)
        {
            const Candidate &cand = cands[i];

            // Away from the key words g's recovered state is exactly X^(R), so
            // only the eight key words are recomputed.
            forro::copy_state(g_state, x_final);
            for (const u8 w : synco::FORRO_LAYOUT)
            {
                const vec pinned_key =
                    simd::bxor(simd::bor(simd::andnot(x0[w], simd::set1(cand.pin.zero[w])),
                                         simd::set1(cand.pin.one[w])),
                               simd::set1(cand.flip[w]));
                key_diff[w] = simd::bxor(x0[w], pinned_key);
                g_state[w]  = forro::ff_inv(keystream[w], pinned_key);
            }

            write_segments(g_state, x_final, key_diff, cand.plan);
            forro::VecBackward::range<DIST_STEP, TOTAL_STEP>(g_state);
            fold_g[i] = masked_value(g_state);
        }
    }

    inline void process_group(const vec *x0, const Candidate *cands, size_t n_cand,
                              vec *acc, vec *scratch)
    {
        vec fold_f;
        run_branch(x0, cands, n_cand, fold_f, scratch);
        for (size_t i = 0; i < n_cand; ++i)
            acc[i] = simd::add(acc[i],
                               simd::bxor(simd::parity(simd::bxor(fold_f, scratch[i])),
                                          simd::set1(1)));
    }

    // Returns median_k |correlation| per candidate. eps'_a is signed per key, so
    // pooling across keys would cancel it; each key is reduced separately.
    inline std::vector<double> estimate(const std::vector<Candidate> &cands, u64 keys, u64 ivs,
                                        u64 seed)
    {
        const u64      groups  = std::max<u64>(1, ivs / LANES);
        const u64      n_iv    = groups * LANES;
        const size_t   n_cand  = cands.size();
        const unsigned hw      = std::thread::hardware_concurrency();
        const size_t   threads = hw > 1 ? hw - 1 : 1;

        std::vector<double> out(n_cand, 0.0);
        if (n_cand == 0 || keys == 0)
            return out;

        std::vector<double>            corr(static_cast<size_t>(keys) * n_cand, 0.0);
        std::vector<std::future<void>> futures;

        for (size_t t = 0; t < threads; ++t)
        {
            const u64 lo = keys * t / threads, hi = keys * (t + 1) / threads;
            if (hi <= lo)
                continue;
            futures.emplace_back(std::async(std::launch::async,
                [&cands, &corr, seed, lo, hi, groups, n_iv, n_cand]()
            {
                std::vector<vec> acc(n_cand), scratch(n_cand);
                std::vector<u64> matches(n_cand);
                vec              x0[STATE_WORDS];
                u32              key[KEY_COUNT];

                for (u64 k = lo; k < hi; ++k)
                {
                    derive_key(seed, k, key);
                    VecRng rng(seed ^ (0xD1B54A32D192ED03ULL * (k + 1)));
                    std::fill(acc.begin(), acc.end(), simd::set1(0));
                    std::fill(matches.begin(), matches.end(), 0ULL);

                    constexpr u64 DRAIN = 256;   // bounds the u32 lane accumulators
                    for (u64 s = 0; s < groups; s += DRAIN)
                    {
                        const u64 count = std::min<u64>(DRAIN, groups - s);
                        for (u64 g = 0; g < count; ++g)
                        {
                            init_state(x0, key, rng);
                            process_group(x0, cands.data(), n_cand, acc.data(), scratch.data());
                        }
                        for (size_t i = 0; i < n_cand; ++i)
                        {
                            matches[i] += simd::hsum(acc[i]);
                            acc[i] = simd::set1(0);
                        }
                    }
                    for (size_t i = 0; i < n_cand; ++i)
                        corr[static_cast<size_t>(k) * n_cand + i] =
                            2.0 * double(matches[i]) / double(n_iv) - 1.0;
                }
            }));
        }
        for (auto &f : futures)
            f.get();

        std::vector<double> column(keys);
        for (size_t i = 0; i < n_cand; ++i)
        {
            for (u64 k = 0; k < keys; ++k)
                column[k] = std::fabs(corr[static_cast<size_t>(k) * n_cand + i]);
            std::nth_element(column.begin(), column.begin() + keys / 2, column.end());
            out[i] = column[keys / 2];
        }
        return out;
    }

    // ---------------- helpers ----------------

    inline synco::Plan plan_for(const std::vector<u16> &set)
    {
        std::array<u32, STATE_WORDS> m{};
        for (u16 idx : set)
        {
            u16 w = 0, b = 0;
            forro::key_bit_to_word(idx, w, b);
            m[w] |= (1u << b);
        }
        return synco::build_plan(m);
    }

    inline std::vector<u16> restricted_key_bits(const synco::Plan &plan,
                                                const std::vector<char> &active)
    {
        std::vector<u16> kr;
        for (size_t c = 0; c < plan.conditions.size(); ++c)
            if (active[c])
                kr.push_back(forro::word_bit_to_key_bit(plan.conditions[c].word,
                                                      plan.conditions[c].restricted_bit));
        std::sort(kr.begin(), kr.end());
        kr.erase(std::unique(kr.begin(), kr.end()), kr.end());
        return kr;
    }

    inline bool read_set_file(const std::string &path, std::vector<u16> &out)
    {
        std::ifstream fin(path);
        if (!fin.is_open())
        {
            std::cerr << "[ERROR] could not open '" << path << "'\n";
            return false;
        }

        std::string body, line;
        while (std::getline(fin, line))
        {
            const size_t hash = line.find('#');
            if (hash != std::string::npos)
                line.resize(hash);
            for (char &ch : line)
                if (ch == ',' || ch == '[' || ch == ']')
                    ch = ' ';
            body += line + ' ';
        }

        std::istringstream iss(body);
        int                b = 0;
        while (iss >> b)
        {
            if (b < 0 || b >= int(NUM_KEY_BITS))
            {
                std::cerr << "[ERROR] key bit " << b << " out of range\n";
                return false;
            }
            out.push_back(u16(b));
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return !out.empty();
    }

    inline std::string list_str(const std::vector<u16> &v)
    {
        std::string s = "[";
        for (size_t i = 0; i < v.size(); ++i)
            s += (i ? ", " : "") + std::to_string(v[i]);
        return s + "]";
    }

    inline std::string mask_str()
    {
        std::string s;
        for (size_t i = 0; i < OUTPUT_MASK_BITS.size(); ++i)
            s += (i ? " + " : "") + std::string("v") + std::to_string(OUTPUT_MASK_BITS[i].first)
                 + "[" + std::to_string(OUTPUT_MASK_BITS[i].second) + "]";
        return s;
    }

    // ---------------- driver ----------------

    // drop: indices of conditions that are inessential and removed. They are
    // printed with a '*' and excluded from theta and k^R.
    inline int run(int argc, char **argv, std::initializer_list<size_t> drop = {})
    {
        std::string path     = FORRO_CPNB_FILE;
        u64         seed     = 1;
        int         keys     = FORRO_DEFAULT_LOG2_KEYS, ivs = FORRO_DEFAULT_LOG2_IVS;

        for (int i = 1; i < argc; ++i)
        {
            const std::string a = argv[i];
            if      (a.rfind("--keys=", 0) == 0) keys = atoi(a.c_str() + 7);
            else if (a.rfind("--ivs=",  0) == 0) ivs  = atoi(a.c_str() + 6);
            else if (a.rfind("--seed=", 0) == 0) seed = strtoull(a.c_str() + 7, nullptr, 10);
            else if (a == "-h" || a == "--help")
            {
                printf("usage: %s [FILE] [--keys=N] [--ivs=N] [--seed=N]\n"
                       "  FILE      the CPNB set (default %s)\n"
                       "  --keys=N  log2 number of random keys   (default %d)\n"
                       "  --ivs=N   log2 number of IVs per key   (default %d)\n"
                       "  --seed=N  PRNG seed                    (default 1)\n",
                       argv[0], FORRO_CPNB_FILE, FORRO_DEFAULT_LOG2_KEYS,
                       FORRO_DEFAULT_LOG2_IVS);
                return 0;
            }
            else if (a.rfind("--", 0) != 0) path = a;
        }
        if (keys < 1 || keys > 30 || ivs < 1 || ivs > 40)
        {
            fprintf(stderr, "[ERROR] --keys and --ivs are log2 counts, out of range\n");
            return 1;
        }

        if (!forro::self_check())
        {
            fprintf(stderr, "[ERROR] the Forro self-check failed, refusing to measure\n");
            return 1;
        }

        std::vector<u16> set;
        if (!read_set_file(path, set))
            return 1;

        const std::time_t t  = std::time(nullptr);
        const std::tm    *lt = std::localtime(&t);
        printf("######## %d-%d-%d %d:%02d:%02d ########\n\n",
               lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min,
               lt->tm_sec);

        const auto t_begin = std::chrono::steady_clock::now();

        printf("Attacking %.2f-round Forro with a %.0f-round differential-linear distinguisher\n"
               "and a %.2f-round backward approximation equipped with the syncopation\n"
               "technique.\n\n",
               TOTAL_ROUNDS, DISTINGUISHING_ROUND, TOTAL_ROUNDS - DISTINGUISHING_ROUND);

        printf(" ###########################  BEGIN  #################################\n\n");

        printf("Forward differential-linear distinguisher of %.0f rounds:\n", DISTINGUISHING_ROUND);
        printf("\tGamma_out = %s, with correlation 2^{%.4f}.\n\n\n", mask_str().c_str(),
               FORWARD_CORR_LOG2);

        printf("Backward approximation with the syncopation technique over %d subrounds:\n",
               TOTAL_STEP - DIST_STEP);
        printf("\tWith %zu CPNBs read from '%s',\n\n", set.size(), path.c_str());
        printf("\t- self-check: %d rounds forward and backward, SIMD == scalar -- OK\n",
               14);
        printf("\t- %s, %zu lanes\n", SIMD_ISA, LANES);
        printf("\t- number of keys: 2^{%d}\n", keys);
        printf("\t- number of IVs:  2^{%d}\n", ivs);
        printf("\t- PRNG seed: %llu\n", (unsigned long long)seed);
        printf("\t- g pins each CPNB block as %s\n\n",
               USE_PATTERN_SPLIT ? "100...0" : "all zero");

        init_masks();

        // Algorithm 1: syncopated segments and restricted segments.
        synco::Plan       plan = plan_for(set);
        std::vector<char> active(plan.conditions.size(), 1);
        for (size_t c : drop)
        {
            if (c >= active.size())
            {
                fprintf(stderr, "[ERROR] condition %zu to drop does not exist\n", c);
                return 1;
            }
            active[c] = 0;
        }
        if (drop.size())
            synco::apply_active(plan, active);

        size_t theta = 0;
        for (char x : active)
            theta += x ? 1 : 0;

        // Algorithm 1's conditions, theta and k^R, off by default as in the
        // reference artifact. Set to 1 to print them.
        constexpr int flagDetails = 0;
        if (flagDetails == 1)
        {
            printf("\t- Algorithm 1 gives %zu conditions (* = inessential, removed):\n",
                   plan.conditions.size());
            for (size_t i = 0; i < plan.conditions.size(); ++i)
            {
                const synco::Condition &c = plan.conditions[i];
                printf("\t    %2u%c x[%2u] bit %2u (key bit %3u) -> seg %08x carry %08x\n",
                       c.index, active[i] ? ':' : '*', c.word, c.restricted_bit,
                       forro::word_bit_to_key_bit(c.word, c.restricted_bit),
                       c.seg_mask, c.carry_mask);
            }
            const std::vector<u16> kr = restricted_key_bits(plan, active);
            printf("\n\t- theta = %zu, |k^R| = %zu\n", theta, kr.size());
            printf("\t- k^R = %s\n\n", list_str(kr).c_str());
        }

        Candidate cand;
        cand.pin  = synco::assign(set, &forro::key_bit_to_word, USE_PATTERN_SPLIT);
        cand.plan = plan;

        /*According to Lemmas 1 and 2, under the condition that equations are
          satisfied, the syncopated segments are independent with PNBs and known.
          Thus, theoretically estimate the backward correlation on the premise
          that the values of syncopated segments are known.*/
        const double eps = estimate({cand}, 1ULL << keys, 1ULL << ivs, seed)[0];
        const double eps_a = eps * eps;

        printf("\t- |eps'_a| = %.6f (2^{%.4f})\n", eps, std::log2(eps));
        printf("\t- eps_a = (eps'_a)^2 = %.8f (2^{%.4f})\n", eps_a, std::log2(eps_a));

        printf("\n #############################  END  ###############################\n\n\n");

        const double secs =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_begin).count();
        printf("Totally used time: %lf seconds.\n\n", secs);
        return 0;
    }
}
