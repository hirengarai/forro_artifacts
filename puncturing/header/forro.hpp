// Forro, 32-bit reference implementation.

#pragma once

#include <bit>
#include <cstdint>

namespace forro {

inline constexpr uint32_t CONSTANT[4] = {0x746C6F76, 0x61616461,
                                         0x72626173, 0x61636E61};
inline constexpr int CONSTANT_AT[4] = {6, 7, 14, 15};
inline constexpr int KEY_AT[8]      = {0, 1, 2, 3, 8, 9, 10, 11};
inline constexpr int IV_AT[4]       = {4, 5, 12, 13};

inline constexpr int COL[4][5] = {
    {0, 4, 8, 12, 3}, {1, 5, 9, 13, 0}, {2, 6, 10, 14, 1}, {3, 7, 11, 15, 2}};
inline constexpr int DIAG[4][5] = {
    {0, 5, 10, 15, 3}, {1, 6, 11, 12, 0}, {2, 7, 8, 13, 1}, {3, 4, 9, 14, 2}};

inline void srf(uint32_t &a, uint32_t &b, uint32_t &c, uint32_t &d,
                uint32_t &e) {
  d += e;  c ^= d;  b = std::rotl(uint32_t(b + c), 10);
  a += b;  e ^= a;  d = std::rotl(uint32_t(d + e), 27);
  c += d;  b ^= c;  a = std::rotl(uint32_t(a + b), 8);
}

inline void init(uint32_t *x, const uint32_t *key, const uint32_t *iv) {
  for (int i = 0; i < 4; ++i) x[CONSTANT_AT[i]] = CONSTANT[i];
  for (int i = 0; i < 8; ++i) x[KEY_AT[i]] = key[i];
  for (int i = 0; i < 4; ++i) x[IV_AT[i]] = iv[i];
}

// Applies subrounds [from, to). One subround is one SRF on one tuple, so
// four of them make a round; odd rounds use COL and even rounds DIAG.
inline void run(uint32_t *x, int from, int to) {
  for (int s = from; s < to; ++s) {
    const int *t = ((s / 4) % 2 == 0) ? COL[s % 4] : DIAG[s % 4];
    srf(x[t[0]], x[t[1]], x[t[2]], x[t[3]], x[t[4]]);
  }
}

}  // namespace forro
