// Shared machinery for the puncturing verification programs: the .punc
// file, the punctured map g, and the reporting helpers.

#pragma once

#include "forro.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace punc {

inline std::vector<uint32_t> without(const std::vector<uint32_t> &bits,
                                     uint32_t mask);
inline uint32_t gather_bits(uint32_t value,
                            const std::vector<uint32_t> &positions);
inline uint32_t random_word();

/// The bits of `value` selected by `mask`, packed down to the low end.
inline uint32_t squeeze(uint32_t value, uint32_t mask) {
  uint32_t out = 0, at = 0;
  while (mask) {
    const uint32_t bit = mask & (~mask + 1u);
    out |= ((value & bit) ? 1u : 0u) << at++;
    mask ^= bit;
  }
  return out;
}
inline void transform(double *row, size_t width);

struct KeyBranch;
inline std::vector<double> build_branch_table(const KeyBranch &branch,
                                              uint32_t key_word_value,
                                              bool key_independent_only);
inline std::vector<double> fold_branch(const KeyBranch &branch,
                                       uint32_t key_word_value,
                                       bool key_independent_only);

// A state word where the key is added: the attacker sees z = v + k and has
// to guess part of k to get at v.
struct KeyBranch {
  uint32_t state_word = 0;   // which of Forro's 16 state words
  uint32_t key_word   = 0;   // which of the 8 key words is added to it
  std::vector<uint32_t> keystream_bits;  // the bits of z the attack reads
  std::vector<uint32_t> key_bits;        // the bits of k the attack guesses
  std::vector<uint32_t> masks;           // this branch's distinct masks

  // The branch spectrum, as parallel arrays: entry i says mask
  // masks[mask_index[i]] contributes coefficient[i] on z_mask[i], k_mask[i].
  std::vector<uint32_t> mask_index, z_mask, k_mask;
  std::vector<double>   coefficient;

  // A bit set in EVERY mask only flips a fixed sign, so it is never
  // guessed. That is why the count is 60 and not 62.
  std::vector<uint32_t> guessed_key_bits;   // key_bits minus the sign bits
  std::vector<uint32_t> sign_key_bits;      // the ones that only flip a sign
  std::vector<uint32_t> varying_keystream_bits;

  size_t table_width() const { return size_t(1) << keystream_bits.size(); }
};

struct PuncturedMap {
  uint64_t n_terms = 0;
  double   rho2_before_key = 0;
  double   rho2_claimed    = 0;
  uint32_t target_word = 0, target_bit = 0;
  uint32_t target_subround = 0, total_subrounds = 0;

  std::vector<KeyBranch> key_branch;
  // No key added, so these words are known outright.
  std::vector<uint32_t>  known_word;
  size_t known_bits_read = 0, known_bits_varying = 0;   // derived

  // The hull, one row per term.
  std::vector<uint32_t> term_mask_index;
  std::vector<uint32_t> term_known_mask;
  std::vector<double>   term_coefficient;
};


/// gather_bits(0b1011, {0, 3}) is 0b11.
inline uint32_t gather_bits(uint32_t value,
                            const std::vector<uint32_t> &positions) {
  uint32_t index = 0;
  for (size_t j = 0; j < positions.size(); ++j)
    index |= ((value >> positions[j]) & 1u) << j;
  return index;
}

/// One branch's contribution for a fixed key, as a table of VALUES, row m
/// holding mask masks[m] indexed by the keystream bits the branch reads.
///
/// Fixing the key turns each key mask into a sign, leaving a function of
/// the kept keystream bits -- but leaving it as a Walsh SPECTRUM. The
/// sampling needs the function, so each row is transformed.
/// `guessed_only` keeps just the entries whose key mask is zero on the
/// bits the attack guesses. Averaging g over all wrong guesses leaves
/// exactly those, so this builds the floor the wrong-key column sits at.
inline std::vector<double> fold_branch(const KeyBranch &branch,
                                       uint32_t key_word_value,
                                       bool key_independent_only = false) {
  uint32_t guessed = 0;
  for (uint32_t bit : branch.guessed_key_bits) guessed |= 1u << bit;

  const size_t width = branch.table_width();
  std::vector<double> spectrum(branch.masks.size() * width, 0.0);
  for (size_t i = 0; i < branch.coefficient.size(); ++i) {
    if (key_independent_only && (branch.k_mask[i] & guessed)) continue;
    const bool negative =
        std::popcount(branch.k_mask[i] & key_word_value) & 1;
    const size_t at = branch.mask_index[i] * width
                      + gather_bits(branch.z_mask[i], branch.keystream_bits);
    spectrum[at] += negative ? -branch.coefficient[i] : branch.coefficient[i];
  }
  return spectrum;
}

/// Walsh-Hadamard over `width` entries, in place: coefficients -> values.
inline void transform(double *row, size_t width) {
  for (size_t span = 1; span < width; span *= 2)
    for (size_t start = 0; start < width; start += 2 * span)
      for (size_t i = start; i < start + span; ++i) {
        const double lo = row[i], hi = row[i + span];
        row[i]        = lo + hi;
        row[i + span] = lo - hi;
      }
}

inline std::vector<double> build_branch_table(const KeyBranch &branch,
                                              uint32_t key_word_value,
                                              bool key_independent_only
                                                  = false) {
  std::vector<double> table =
      fold_branch(branch, key_word_value, key_independent_only);
  const size_t width = branch.table_width();
  for (size_t m = 0; m < branch.masks.size(); ++m)
    transform(table.data() + m * width, width);
  return table;
}

/// A sample's randomness is a function of (seed, key, sample) alone, so a
/// run reproduces whatever the thread count is and however the work splits.
inline uint64_t mix(uint64_t z) {
  z += 0x9E3779B97F4A7C15ULL;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

inline void draw_words(uint64_t seed, uint32_t *out, size_t count) {
  for (size_t i = 0; i < count; ++i) {
    seed = mix(seed);
    out[i] = uint32_t(seed >> 32);
  }
}

inline uint64_t sample_seed(uint64_t seed, uint32_t key_index,
                            uint64_t sample) {
  return mix(mix(seed ^ (uint64_t(key_index) << 40)) ^ sample);
}

inline thread_local uint64_t rng_state = std::random_device{}();

inline void seed_stream(uint64_t seed) { rng_state = seed; }

inline std::vector<uint32_t> without(const std::vector<uint32_t> &bits,
                                     uint32_t mask) {
  std::vector<uint32_t> kept;
  for (uint32_t bit : bits)
    if (!((mask >> bit) & 1u)) kept.push_back(bit);
  return kept;
}

inline uint32_t random_word() {
  uint64_t z = (rng_state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return uint32_t(z ^ (z >> 31));
}


template <typename T> inline T read_value(std::istream &in) {
  T value{};
  in.read(reinterpret_cast<char *>(&value), sizeof(T));
  return value;
}

template <typename T>
inline std::vector<T> read_array(std::istream &in, size_t count) {
  std::vector<T> values(count);
  if (count)
    in.read(reinterpret_cast<char *>(values.data()),
            std::streamsize(count * sizeof(T)));
  return values;
}

inline PuncturedMap read_map(const char *path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    std::cerr << "cannot open " << path << "\n";
    std::exit(1);
  }
  char magic[8];
  in.read(magic, 8);
  if (std::memcmp(magic, "FORROPUN", 8) != 0
      || read_value<uint32_t>(in) != 1u) {
    std::cerr << path << " is not a version 1 .punc file\n";
    std::exit(1);
  }
  if (read_value<uint32_t>(in) != 32u) {
    std::cerr << path << " is not for 32-bit words\n";
    std::exit(1);
  }

  PuncturedMap map;
  const uint32_t n_key_branches = read_value<uint32_t>(in);
  const uint32_t n_known_words  = read_value<uint32_t>(in);
  if (n_key_branches == 0 || n_key_branches > 8 || n_known_words > 8) {
    std::cerr << path << " declares " << n_key_branches << " key and "
              << n_known_words << " known branches; 1..8 and 0..8 are the "
              << "limits this build holds\n";
    std::exit(1);
  }
  map.n_terms         = read_value<uint64_t>(in);
  map.rho2_before_key = read_value<double>(in);
  map.rho2_claimed    = read_value<double>(in);
  map.target_word     = read_value<uint32_t>(in);
  map.target_bit      = read_value<uint32_t>(in);
  map.target_subround = read_value<uint32_t>(in);
  map.total_subrounds = read_value<uint32_t>(in);

  map.key_branch.resize(n_key_branches);
  for (KeyBranch &branch : map.key_branch) {
    branch.state_word = read_value<uint32_t>(in);
    branch.key_word   = read_value<uint32_t>(in);
    branch.keystream_bits =
        read_array<uint32_t>(in, read_value<uint32_t>(in));
    branch.key_bits = read_array<uint32_t>(in, read_value<uint32_t>(in));
    branch.masks    = read_array<uint32_t>(in, read_value<uint32_t>(in));
    const uint64_t n = read_value<uint64_t>(in);
    branch.mask_index  = read_array<uint32_t>(in, n);
    branch.z_mask      = read_array<uint32_t>(in, n);
    branch.k_mask      = read_array<uint32_t>(in, n);
    branch.coefficient = read_array<double>(in, n);
  }

  for (KeyBranch &branch : map.key_branch) {
    if (branch.masks.empty() || branch.keystream_bits.size() > 24
        || branch.state_word > 15 || branch.key_word > 7) {
      std::cerr << path << " has an out-of-range branch\n";
      std::exit(1);
    }
    uint32_t always = 0xFFFFFFFFu;
    for (uint32_t mask : branch.masks) always &= mask;
    branch.guessed_key_bits       = without(branch.key_bits, always);
    branch.sign_key_bits          = without(branch.key_bits,
                                            ~always);
    branch.varying_keystream_bits = without(branch.keystream_bits, always);
  }

  map.known_word = read_array<uint32_t>(in, n_known_words);
  map.term_mask_index =
      read_array<uint32_t>(in, size_t(map.n_terms) * n_key_branches);
  map.term_known_mask =
      read_array<uint32_t>(in, size_t(map.n_terms) * n_known_words);
  map.term_coefficient = read_array<double>(in, size_t(map.n_terms));

  if (!in) {
    std::cerr << path << " is truncated\n";
    std::exit(1);
  }

  if (map.target_word > 15 || map.target_bit > 31
      || map.target_subround > map.total_subrounds
      || map.total_subrounds == 0 || map.total_subrounds > 256
      || map.n_terms == 0) {
    std::cerr << path << " has an out-of-range header\n";
    std::exit(1);
  }
  for (uint64_t i = 0; i < uint64_t(map.n_terms) * n_key_branches; ++i)
    if (map.term_mask_index[i]
        >= map.key_branch[i % n_key_branches].masks.size()) {
      std::cerr << path << " has a term pointing outside its branch\n";
      std::exit(1);
    }
  // A word the map calls known must not be one the key is added to, or the
  // recovery keystream - initial would be using a secret.
  for (uint32_t w : map.known_word) {
    if (w > 15) {
      std::cerr << path << " names state word " << w << "\n";
      std::exit(1);
    }
    for (int k = 0; k < 8; ++k)
      if (forro::KEY_AT[k] == int(w)) {
        std::cerr << path << " calls v" << w
                  << " known, but the key is added to it\n";
        std::exit(1);
      }
  }

  // What g reads of a known word is just its active bits.
  for (size_t j = 0; j < n_known_words; ++j) {
    uint32_t any = 0, all = 0xFFFFFFFFu;
    for (uint64_t term = 0; term < map.n_terms; ++term) {
      const uint32_t mask = map.term_known_mask[term * n_known_words + j];
      any |= mask;
      all &= mask;
    }
    map.known_bits_read    += std::popcount(any);
    map.known_bits_varying += std::popcount(uint32_t(any & ~all));
  }
  return map;
}


inline std::string as_power_of_two(double value) {
  if (value == 0) return "0";
  char buffer[32];
  // Two decimals, as the paper quotes them.
  std::snprintf(buffer, sizeof(buffer), "%s2^%.2f", value < 0 ? "-" : "",
                std::log2(std::fabs(value)));
  return buffer;
}

/// printf pads by bytes; the labels hold multi-byte characters.
inline std::string padded(const std::string &text, size_t columns) {
  size_t width = 0;
  for (unsigned char c : text)
    if ((c & 0xC0) != 0x80) ++width;
  return text + std::string(columns > width ? columns - width : 0, ' ');
}

inline void rule() {
  std::printf(" ");
  for (int i = 0; i < 62; ++i) std::printf("\u2500");
  std::printf("\n");
}

inline void print_line(const char *label, const std::string &value,
                       const char *note = "") {
  std::printf(" %s", padded(label, 18).c_str());
  if (*note)
    std::printf("%s  %s", padded(value, 16).c_str(), note);
  else
    std::printf("%s", value.c_str());
  std::printf("\n");
}

inline std::string grouped(uint64_t value) {
  std::string digits = std::to_string(value), out;
  for (size_t i = 0; i < digits.size(); ++i) {
    if (i && (digits.size() - i) % 3 == 0) out += ',';
    out += digits[i];
  }
  return out;
}


// Evaluates g for one key. Each thread holds its own, so the small branch
// tables and the per-sample scratch never cross threads.
class Evaluator {
 public:
  /// With `key_independent_only`, g is averaged over every wrong guess --
  /// which is the floor the wrong-key control actually sits at, and not
  /// zero as one might assume.
  explicit Evaluator(const PuncturedMap &map, const uint32_t *key,
                     bool key_independent_only = false)
      : map_(map), table_(map.key_branch.size()),
        column_(map.key_branch.size()) {
    for (size_t b = 0; b < map.key_branch.size(); ++b) {
      const KeyBranch &branch = map.key_branch[b];
      table_[b] = build_branch_table(branch, key[branch.key_word],
                                     key_independent_only);
      column_[b].resize(branch.masks.size());
    }
  }

  /// `initial` is the public part of the starting state -- constants and
  /// IV -- and `keystream` what the encryption produced. A known word's
  /// value is recovered as keystream - initial, exactly as an attacker
  /// would; nothing here ever touches the cipher's internal state.
  double operator()(const uint32_t *initial, const uint32_t *keystream) {
    const size_t n_key = map_.key_branch.size();
    const size_t n_known = map_.known_word.size();

    for (size_t b = 0; b < n_key; ++b) {
      const KeyBranch &branch = map_.key_branch[b];
      const uint32_t at = gather_bits(keystream[branch.state_word],
                                      branch.keystream_bits);
      const size_t width = branch.table_width();
      for (size_t m = 0; m < branch.masks.size(); ++m)
        column_[b][m] = table_[b][m * width + at];
    }

    uint32_t known[8];
    for (size_t j = 0; j < n_known; ++j) {
      const uint32_t w = map_.known_word[j];
      known[j] = keystream[w] - initial[w];
    }

    double g = 0;
    for (uint64_t term = 0; term < map_.n_terms; ++term) {
      const uint32_t *mask_index =
          map_.term_mask_index.data() + term * n_key;
      const uint32_t *known_mask =
          map_.term_known_mask.data() + term * n_known;
      uint32_t bits = 0;
      for (size_t j = 0; j < n_known; ++j) bits ^= known_mask[j] & known[j];
      double value = (std::popcount(bits) & 1)
                         ? -map_.term_coefficient[term]
                         : map_.term_coefficient[term];
      for (size_t b = 0; b < n_key; ++b) value *= column_[b][mask_index[b]];
      g += value;
    }
    return g;
  }

 private:
  const PuncturedMap &map_;
  std::vector<std::vector<double>> table_, column_;
};

/// A wrong guess differs from the truth only in bits the attack guesses.
inline void make_wrong_key(const PuncturedMap &map, const uint32_t *key,
                           uint32_t *wrong) {
  std::memcpy(wrong, key, 8 * sizeof(uint32_t));
  for (const KeyBranch &branch : map.key_branch)
    for (uint32_t bit : branch.guessed_key_bits)
      if (random_word() & 1) wrong[branch.key_word] ^= 1u << bit;
  if (std::memcmp(key, wrong, 8 * sizeof(uint32_t)) == 0)
    wrong[map.key_branch[0].key_word] ^=
        1u << map.key_branch[0].guessed_key_bits.front();
}

/// The mean of `values`, and the standard error of that mean.
inline void summarise(const std::vector<double> &values, double &mean,
                      double &error) {
  mean = 0;
  for (double v : values) mean += v;
  mean /= double(values.size());
  double spread = 0;
  for (double v : values) spread += (v - mean) * (v - mean);
  error = values.size() > 1
              ? std::sqrt(spread / double(values.size() - 1)
                          / double(values.size()))
              : 0.0;
}

/// The block every experiment prints before it starts.
inline void print_map(const PuncturedMap &map, const char *path,
                      const char *sampling) {
  size_t guessed = 0, sign_bits = 0, read = 0, varying = 0;
  for (const KeyBranch &branch : map.key_branch) {
    guessed   += branch.guessed_key_bits.size();
    sign_bits += branch.sign_key_bits.size();
    read      += branch.keystream_bits.size();
    varying   += branch.varying_keystream_bits.size();
  }
  read    += map.known_bits_read;
  varying += map.known_bits_varying;

  char target[96], attack[64], reads[96], guesses[96];
  std::snprintf(target, sizeof(target), "v%u[%u] after subround %u",
                map.target_word, map.target_bit, map.target_subround);
  std::snprintf(attack, sizeof(attack), "%u subrounds (%.2f rounds)",
                map.total_subrounds, 0.25 * map.total_subrounds);
  std::snprintf(reads, sizeof(reads), "%zu keystream bits (%zu varying)",
                read, varying);
  std::snprintf(guesses, sizeof(guesses),
                "%zu key bits (%zu more only flip the sign)", guessed,
                sign_bits);

  print_line("map", path);
  print_line("target", target);
  print_line("attack", attack);
  print_line("hull", grouped(map.n_terms) + " terms");
  print_line("g reads", reads);
  print_line("g guesses", guesses);
  print_line("sampling", sampling);
}

}  // namespace punc
