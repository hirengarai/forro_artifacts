// Experiment 2: does the puncturing loss on a pair factorise?
//
// This is "Verification of Puncturing Correlation, rho^4" from Appendix C.1
// of Florez-Gutierrez and Todo, "Improved Cryptanalysis of ChaCha: Beating
// PNBs with Bit Puncturing", EUROCRYPT 2025,
// https://eprint.iacr.org/2025/437 , carried over to Forro.
//
// A differential-linear attack analyses two texts and prices the pair at
// rho^4, which assumes the two are independent. They share a key and are
// related by the input difference, so that needs checking. Each key is
// therefore compared against itself: one pass measures the single-text
// correlation c on each member and the pair correlation p, and
// independence says p = c * c'. Averaging rho^2 over keys first would
// bias the comparison, so no key average enters the test.
//
//   g++ -std=c++20 -O3 -march=native -pthread -o rho4 rho4.cpp
//   ./rho4 data/f_1_5round.punc [seed]
//   -v anywhere in the arguments adds the header, controls and verdict

#include "header/punc.hpp"

#include <chrono>
#include <thread>

// --------- Control Panel ---------

static constexpr int KEYS            = 1024;
static constexpr int SAMPLES_PER_KEY = 1 << 12;

// The attack's own pairs: the second text carries the input difference of
// the distinguisher. Set false to use two unrelated texts instead, which
// tests only the shared key and not the relation.
static constexpr bool   ATTACK_PAIRS = true;
static constexpr uint32_t DIFF_WORD  = 5, DIFF_BIT = 18;

static constexpr unsigned THREADS = 0;        // 0 = as many as there are cores
static constexpr uint64_t DEFAULT_SEED = 0;   // 0 = a fresh seed each run
static constexpr const char *DEFAULT_MAP = "data/f_1_5round.punc";

// --------- Types ---------

struct Totals {
  double   first = 0, second = 0, pair = 0;   // f g,  f' g',  f g f' g'
  double   wrong_pair = 0;
  uint64_t count = 0;
};

// --------- Declarations ---------

static Totals take_samples(const punc::PuncturedMap &map, const uint32_t *key,
                           const uint32_t *wrong_key, uint64_t first,
                           uint64_t last, uint32_t key_index, uint64_t seed);
static void one_text(const punc::PuncturedMap &map, const uint32_t *key,
                     const uint32_t *iv, bool flip,
                     punc::Evaluator &right, punc::Evaluator &wrong,
                     double &f, double &g, double &g_wrong);

// --------- Main ---------

// Off by default: the run prints the two measurements, and -v adds the
// header, the per-key table, the control and the verdict.
static bool verbose = false;

int main(int argc, char **argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  std::vector<const char *> args;
  for (int i = 1; i < argc; ++i) {
    if (!std::strcmp(argv[i], "-v")) verbose = true;
    else args.push_back(argv[i]);
  }
  const char *path = args.size() > 0 ? args[0] : DEFAULT_MAP;
  uint64_t seed = args.size() > 1 ? std::strtoull(args[1], nullptr, 0)
                                  : DEFAULT_SEED;
  if (seed == 0)
    seed = (uint64_t(std::random_device{}()) << 32) ^ std::random_device{}();

  const punc::PuncturedMap map = punc::read_map(path);
  unsigned threads = THREADS ? THREADS : std::thread::hardware_concurrency();
  if (threads < 1) threads = 1;

  const auto started = std::chrono::steady_clock::now();
  std::vector<double> gap_per_key, ratio_per_key, wrong_per_key, p_per_key;
  double sum_c = 0, sum_p = 0, sum_cc = 0;

  const bool list_keys = verbose && KEYS <= 20;
  if (list_keys)
    std::printf(" %5s  %12s  %12s  %12s\n", "key", "c", "p", "p - c*c'");
  for (int k = 0; k < KEYS; ++k) {
    uint32_t key[8], wrong_key[8];
    punc::draw_words(punc::mix(seed ^ (uint64_t(k) << 20)), key, 8);
    punc::seed_stream(punc::mix(seed ^ (uint64_t(k) << 21)));
    punc::make_wrong_key(map, key, wrong_key);

    std::vector<Totals> parts(threads);
    std::vector<std::thread> workers;
    for (unsigned t = 0; t < threads; ++t)
      workers.emplace_back([&, t] {
        const uint64_t lo = uint64_t(SAMPLES_PER_KEY) * t / threads;
        const uint64_t hi = uint64_t(SAMPLES_PER_KEY) * (t + 1) / threads;
        parts[t] = take_samples(map, key, wrong_key, lo, hi, uint32_t(k),
                                seed);
      });
    for (std::thread &worker : workers) worker.join();

    Totals sum;
    for (const Totals &part : parts) {
      sum.first  += part.first;   sum.second += part.second;
      sum.pair   += part.pair;    sum.wrong_pair += part.wrong_pair;
      sum.count  += part.count;
    }
    const double n = double(sum.count);
    const double c = sum.first / n, c2 = sum.second / n, p = sum.pair / n;
    gap_per_key.push_back(p - c * c2);
    ratio_per_key.push_back(c * c2 != 0 ? p / (c * c2) : 0.0);
    wrong_per_key.push_back(sum.wrong_pair / n);
    p_per_key.push_back(p);
    sum_c += 0.5 * (c + c2);
    sum_p += p;
    sum_cc += c * c2;
    if (list_keys)
      std::printf(" %5d  %12.6f  %12.6f  %12.6f\n", k, c, p,
                  gap_per_key.back());
    else if (verbose && (k + 1) % (KEYS / 10 ? KEYS / 10 : 1) == 0)
      std::printf(" %5d of %d\n", k + 1, KEYS);
  }

  double gap, gap_error, ratio, ratio_error, wrong, wrong_error;
  punc::summarise(gap_per_key, gap, gap_error);
  punc::summarise(ratio_per_key, ratio, ratio_error);
  punc::summarise(wrong_per_key, wrong, wrong_error);

  char gap_note[96], ratio_note[96];
  std::snprintf(gap_note, sizeof(gap_note), "+-%.6f    %+.2f sigma from 0",
                gap_error, gap_error ? gap / gap_error : 0.0);
  std::snprintf(ratio_note, sizeof(ratio_note), "ratio %.3f    %+.2f sigma    +-%.4f"
                " over %d keys", ratio,
                ratio_error ? (ratio - 1.0) / ratio_error : 0.0, ratio_error,
                KEYS);
  char elapsed[32], elapsed_note[80];
  const double secs = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  if (secs >= 3600)     std::snprintf(elapsed, sizeof(elapsed), "%.1f h", secs / 3600);
  else if (secs >= 60)  std::snprintf(elapsed, sizeof(elapsed), "%.1f min", secs / 60);
  else                  std::snprintf(elapsed, sizeof(elapsed), "%.1f s", secs);
  std::snprintf(elapsed_note, sizeof(elapsed_note),
                "%.0f s, %d keys x %d samples on %u threads", secs, KEYS,
                SAMPLES_PER_KEY, threads);

  if (verbose) punc::rule();
  punc::print_line("\u03c1\u2074 predicted",
                   punc::as_power_of_two(sum_cc / KEYS),
                   "c · c', measured per key and averaged");
  punc::print_line("\u03c1\u2074 measured", punc::as_power_of_two(sum_p / KEYS),
                   "the pair correlation");
  punc::print_line("elapsed", elapsed, elapsed_note);
  if (verbose) {
    char mean_c[32], mean_p[32], naive[96];
    std::snprintf(mean_c, sizeof(mean_c), "%.6f", sum_c / KEYS);
    std::snprintf(mean_p, sizeof(mean_p), "%.6f", sum_p / KEYS);
    std::snprintf(naive, sizeof(naive),
                  "%+.1f%% from the prediction; E[c\u00b2] \u2260 (E[c])\u00b2,"
                  " which is why the test is per key",
                  100.0 * (map.rho2_claimed * map.rho2_claimed
                           - sum_cc / KEYS) / (sum_cc / KEYS));
    punc::print_line("agreement", ratio_note);
    double p_mean, p_error;
    punc::summarise(p_per_key, p_mean, p_error);
    const double var_expected = map.rho2_claimed * map.rho2_claimed
                                / double(SAMPLES_PER_KEY);
    const double var_observed = p_error * p_error * double(KEYS);
    char var_note[96];
    std::snprintf(var_note, sizeof(var_note), "across %d keys, %.1fx predicted",
                  KEYS, var_expected ? var_observed / var_expected : 0.0);
    punc::print_line("variance expected", punc::as_power_of_two(var_expected),
                     "\u03c1\u2074 / samples per key");
    punc::print_line("variance observed", punc::as_power_of_two(var_observed),
                     var_note);
    punc::print_line("mean c", mean_c, "per key, the single-text correlation");
    punc::print_line("mean p", mean_p, "per key, the pair");
    punc::print_line("p - c · c'", std::string(1, ' ') + std::to_string(gap),
                     gap_note);
    punc::print_line("(\u03c1\u00b2)\u00b2",
                     punc::as_power_of_two(map.rho2_claimed
                                           * map.rho2_claimed), naive);
  }

  // Only the factorisation gates this. The wrong-key column here is a
  // product of two pair statistics, whose expectation is not zero and is
  // not what this experiment tests; it is reported, not judged.
  const bool agrees = std::fabs(gap) < 3 * gap_error;
  if (verbose) {
    punc::print_line("wrong key", punc::as_power_of_two(wrong),
                     "reported, not a pass/fail");
    punc::rule();
    punc::print_line("result", agrees ? "factorises" : "DOES NOT FACTORISE");
    punc::rule();
  }
  return agrees ? 0 : 1;
}

// --------- Sampling ---------

static void one_text(const punc::PuncturedMap &map, const uint32_t *key,
                     const uint32_t *iv, bool flip,
                     punc::Evaluator &right, punc::Evaluator &wrong,
                     double &f, double &g, double &g_wrong) {
  uint32_t start[16], state[16], keystream[16], public_state[16];
  forro::init(start, key, iv);
  if (flip) start[DIFF_WORD] ^= 1u << DIFF_BIT;
  std::memcpy(public_state, start, sizeof(start));
  for (int i = 0; i < 8; ++i) public_state[forro::KEY_AT[i]] = 0;
  std::memcpy(state, start, sizeof(start));

  forro::run(state, 0, int(map.target_subround));
  f = (state[map.target_word] >> map.target_bit) & 1 ? -1.0 : 1.0;

  forro::run(state, int(map.target_subround), int(map.total_subrounds));
  for (int w = 0; w < 16; ++w) keystream[w] = start[w] + state[w];

  g       = right(public_state, keystream);
  g_wrong = wrong(public_state, keystream);
}

static Totals take_samples(const punc::PuncturedMap &map, const uint32_t *key,
                           const uint32_t *wrong_key, uint64_t first,
                           uint64_t last, uint32_t key_index, uint64_t seed) {
  punc::Evaluator with_right(map, key), with_wrong(map, wrong_key);

  Totals totals;
  uint32_t iv[4], other_iv[4];

  for (uint64_t sample = first; sample < last; ++sample) {
    const uint64_t s = punc::sample_seed(seed, key_index, sample);
    punc::draw_words(s, iv, 4);
    // The attack's pair shares the IV and differs by the input difference;
    // the unrelated variant draws a second IV instead.
    if (ATTACK_PAIRS) std::memcpy(other_iv, iv, sizeof(iv));
    else punc::draw_words(punc::mix(s ^ 0xA5A5A5A5u), other_iv, 4);

    double f1, g1, w1, f2, g2, w2;
    one_text(map, key, iv, false, with_right, with_wrong, f1, g1, w1);
    one_text(map, key, other_iv, ATTACK_PAIRS, with_right, with_wrong, f2, g2,
             w2);

    totals.first      += f1 * g1;
    totals.second     += f2 * g2;
    totals.pair       += (f1 * g1) * (f2 * g2);
    totals.wrong_pair += (f1 * w1) * (f2 * w2);
    ++totals.count;
  }
  return totals;
}
