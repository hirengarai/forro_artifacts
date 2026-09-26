/*
 * forro675_verification.cpp
 *
 * Verification experiment for the conditional backward correlation of the
 * 6.75-round key-recovery attack on Forro.
 *
 * Reproduces:  n = 80 CPNBs, theta = 14, |k^R| = 14, |eps'_a| = 0.00841.
 *
 * by Hiren
 *
 * Compile: g++ -std=c++20 -O3 -march=native -pthread ./forro675_verification.cpp -o run675
 * Run:     ./run675
 */

#define FORRO_TOTAL_ROUNDS_X4     27                        // 6.75 rounds
#define FORRO_CPNB_FILE           "cpnbs_size80_r675.txt"
#define FORRO_DEFAULT_LOG2_KEYS   10
#define FORRO_DEFAULT_LOG2_IVS    24

#include "verify.hpp"

// Algorithm 1 returns 16 conditions here. Conditions 2 and 8 are inessential:
// dropping them leaves theta = 14 and lowers the attack time, because the
// correlation they buy costs more in data (2^-2 per condition) than it returns.
// See Section IX-B.
int main(int argc, char **argv)
{
    return verify::run(argc, argv, {2, 8});
}
