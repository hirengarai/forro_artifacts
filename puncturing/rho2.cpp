// Experiment 1: the puncturing correlation itself.
//
// This is "Verification of the Puncturing Correlation rho^2" from
// Appendix C.1 of Florez-Gutierrez and Todo, "Improved Cryptanalysis of
// ChaCha: Beating PNBs with Bit Puncturing", EUROCRYPT 2025,
// https://eprint.iacr.org/2025/437 , carried over to Forro.
//
// f is the true target bit, read from inside the encryption; g is the
// punctured map, computed from the keystream and a guessed key. Puncturing
// says rho^2 = <f, g>, so the average of f * g has to reproduce the rho^2
// the search reported. The wrong-key column is the control.
//
//   g++ -std=c++20 -O3 -march=native -pthread -o rho2 rho2.cpp
//   ./rho2 data/f_1_5round.punc [seed]
//   -v anywhere in the arguments adds the header, controls and verdict

#include "header/punc.hpp"

#include <chrono>
#include <thread>

// --------- Control Panel ---------

// rho^2 averages over keys as well as inputs, so accuracy comes from KEYS,
// not from the samples inside one key.
static constexpr int KEYS            = 1024;
static constexpr int SAMPLES_PER_KEY = 1 << 12;

static constexpr unsigned THREADS = 0;        // 0 = as many as there are cores
static constexpr uint64_t DEFAULT_SEED = 0;   // 0 = a fresh seed each run
static constexpr const char *DEFAULT_MAP = "data/f_1_5round.punc";

// --------- Types ---------

struct Totals {
  double   sum_right = 0, sum_wrong = 0, sum_g_squared = 0;
  double   sum_floor = 0;   // g averaged over every wrong guess
  uint64_t count = 0;
};

// --------- Declarations ---------

static Totals take_samples(const punc::PuncturedMap &map, const uint32_t *key,
                           const uint32_t *wrong_key, uint64_t first,
                           uint64_t last, uint32_t key_index, uint64_t seed);

// --------- Main ---------

// Off by default: the run prints the claim and the measurement, and -v
// adds the header, the controls and the verdict.
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
  punc::seed_stream(seed);

  const punc::PuncturedMap map = punc::read_map(path);
  unsigned threads = THREADS ? THREADS : std::thread::hardware_concurrency();
  if (threads < 1) threads = 1;
  char sampling[96];
  std::snprintf(sampling, sizeof(sampling),
                "%d keys x %d samples, %u threads", KEYS, SAMPLES_PER_KEY,
                threads);


  const auto started = std::chrono::steady_clock::now();
  const bool list_keys = verbose && KEYS <= 20;
  if (list_keys) std::printf(" %5s  %12s  %12s\n", "key", "right", "wrong");

  std::vector<double> right_per_key, wrong_per_key, g2_per_key,
      floor_per_key;

  for (int k = 0; k < KEYS; ++k) {
    uint32_t key[8], wrong_key[8];
    punc::draw_words(punc::mix(seed ^ (uint64_t(k) << 20)), key, 8);
    punc::seed_stream(punc::mix(seed ^ (uint64_t(k) << 21)));
    punc::make_wrong_key(map, key, wrong_key);

    // Threads take disjoint ranges of the SAME sample list, so the result
    // does not depend on how many there are.
    std::vector<Totals> parts(threads);
    std::vector<std::thread> workers;
    for (unsigned t = 0; t < threads; ++t)
      workers.emplace_back([&, t] {
        const uint64_t first = uint64_t(SAMPLES_PER_KEY) * t / threads;
        const uint64_t last  = uint64_t(SAMPLES_PER_KEY) * (t + 1) / threads;
        parts[t] = take_samples(map, key, wrong_key, first, last,
                                uint32_t(k), seed);
      });
    for (std::thread &worker : workers) worker.join();

    Totals sum;
    for (const Totals &part : parts) {
      sum.sum_right     += part.sum_right;
      sum.sum_wrong     += part.sum_wrong;
      sum.sum_g_squared += part.sum_g_squared;
      sum.sum_floor     += part.sum_floor;
      sum.count         += part.count;
    }
    right_per_key.push_back(sum.sum_right / double(sum.count));
    wrong_per_key.push_back(sum.sum_wrong / double(sum.count));
    g2_per_key.push_back(sum.sum_g_squared / double(sum.count));
    floor_per_key.push_back(sum.sum_floor / double(sum.count));

    if (list_keys)
      std::printf(" %5d  %12.6f  %12.6f\n", k, right_per_key.back(),
                  wrong_per_key.back());
    else if (verbose && (k + 1) % (KEYS / 10 ? KEYS / 10 : 1) == 0)
      std::printf(" %5d of %d\n", k + 1, KEYS);
  }

  double measured, error, wrong, wrong_error, g2, g2_error, floor,
      floor_error;
  punc::summarise(right_per_key, measured, error);
  punc::summarise(wrong_per_key, wrong, wrong_error);
  punc::summarise(g2_per_key, g2, g2_error);
  punc::summarise(floor_per_key, floor, floor_error);

  char ratio[96], wrong_note[64];
  std::snprintf(ratio, sizeof(ratio), "ratio %.3f    %+.2f sigma    +-%.4f"
                " over %d keys", measured / map.rho2_claimed,
                error ? (measured - map.rho2_claimed) / error : 0.0,
                error, KEYS);
  const double wrong_gap = wrong - floor;
  const double wrong_band = std::sqrt(wrong_error * wrong_error
                                      + floor_error * floor_error);
  std::snprintf(wrong_note, sizeof(wrong_note),
                "+-%.4f    %+.2f sigma from the floor", wrong_error,
                wrong_band ? wrong_gap / wrong_band : 0.0);
  char elapsed[32], elapsed_note[80];
  const double secs = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  if (secs >= 3600)     std::snprintf(elapsed, sizeof(elapsed), "%.1f h", secs / 3600);
  else if (secs >= 60)  std::snprintf(elapsed, sizeof(elapsed), "%.1f min", secs / 60);
  else                  std::snprintf(elapsed, sizeof(elapsed), "%.1f s", secs);
  std::snprintf(elapsed_note, sizeof(elapsed_note),
                "%.0f s, %d keys x %d samples on %u threads", secs, KEYS,
                SAMPLES_PER_KEY, threads);

  const bool agrees = std::fabs(measured - map.rho2_claimed) < 3 * error
                      && std::fabs(g2 - map.rho2_claimed) < 3 * g2_error
                      && std::fabs(wrong_gap) < 3 * wrong_band;

  if (verbose) punc::rule();
  punc::print_line("\u03c1\u00b2 claimed",
                   punc::as_power_of_two(map.rho2_claimed));
  punc::print_line("\u03c1\u00b2 measured", punc::as_power_of_two(measured));
  punc::print_line("\u2016g\u2016\u00b2", punc::as_power_of_two(g2),
                   "a self-check on g alone; must equal \u03c1\u00b2");
  punc::print_line("elapsed", elapsed, elapsed_note);
  if (verbose) {
    punc::print_line("agreement", ratio);
    // Their Appendix C.1 reports the variance against ||g||^2 / |S|, the
    // spread one key would show. Ours is the spread ACROSS keys, so an
    // excess is the correlation's key dependence, not sampling noise.
    const double var_expected = g2 / double(SAMPLES_PER_KEY);
    const double var_observed = error * error * double(KEYS);
    char var_note[96];
    std::snprintf(var_note, sizeof(var_note), "across %d keys, %.1fx predicted",
                  KEYS, var_expected ? var_observed / var_expected : 0.0);
    punc::print_line("variance expected", punc::as_power_of_two(var_expected),
                     "\u2016g\u2016\u00b2 / samples per key");
    punc::print_line("variance observed", punc::as_power_of_two(var_observed),
                     var_note);
    punc::print_line("wrong key", punc::as_power_of_two(wrong), wrong_note);
    punc::print_line("wrong-key floor", punc::as_power_of_two(floor),
                     "g averaged over every wrong guess");
    punc::rule();
    punc::print_line("result", agrees ? "agrees" : "DISAGREES");
    punc::rule();
  }

  return agrees ? 0 : 1;
}

// --------- Sampling ---------

static Totals take_samples(const punc::PuncturedMap &map, const uint32_t *key,
                           const uint32_t *wrong_key, uint64_t first,
                           uint64_t last, uint32_t key_index, uint64_t seed) {
  punc::Evaluator with_right(map, key), with_wrong(map, wrong_key);
  // A wrong guess flips a random half of the guessed bits, so terms whose
  // key mask misses them all survive untouched. That leaves a floor, not
  // zero, and this evaluator computes it.
  punc::Evaluator floor(map, key, true);

  Totals totals;
  uint32_t start[16], state[16], keystream[16], public_state[16], iv[4];

  for (uint64_t sample = first; sample < last; ++sample) {
    punc::draw_words(punc::sample_seed(seed, key_index, sample), iv, 4);
    forro::init(start, key, iv);
    // What an attacker sees of the starting state: constants and IV, with
    // the key words left at zero.
    std::memcpy(public_state, start, sizeof(start));
    for (int i = 0; i < 8; ++i) public_state[forro::KEY_AT[i]] = 0;
    std::memcpy(state, start, sizeof(start));

    forro::run(state, 0, int(map.target_subround));
    const double f =
        (state[map.target_word] >> map.target_bit) & 1 ? -1.0 : 1.0;

    forro::run(state, int(map.target_subround), int(map.total_subrounds));
    for (int w = 0; w < 16; ++w) keystream[w] = start[w] + state[w];

    const double g = with_right(public_state, keystream);
    totals.sum_right     += f * g;
    totals.sum_wrong     += f * with_wrong(public_state, keystream);
    totals.sum_floor     += f * floor(public_state, keystream);
    totals.sum_g_squared += g * g;
    ++totals.count;
  }
  return totals;
}
