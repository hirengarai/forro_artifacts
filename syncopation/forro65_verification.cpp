/*
 * forro65_verification.cpp
 *
 * Verification experiment for the conditional backward correlation of the
 * 6.5-round key-recovery attack on Forro.
 *
 * Reproduces:  n = 113 CPNBs, theta = 15, |k^R| = 15, |eps'_a| = 0.0422.
 *
 * by Hiren
 *
 * Compile: g++ -std=c++20 -O3 -march=native -pthread ./forro65_verification.cpp -o run65
 * Run:     ./run65
 */

#define FORRO_TOTAL_ROUNDS_X4     26                        // 6.5 rounds
#define FORRO_CPNB_FILE           "cpnbs_size113_r65.txt"
#define FORRO_DEFAULT_LOG2_KEYS   12
#define FORRO_DEFAULT_LOG2_IVS    18

#include "verify.hpp"

// Algorithm 1 returns 16 conditions here. Condition 9 is inessential: dropping it
// leaves theta = 15 and lowers the attack time, because the correlation it buys
// costs more in data (2^-2 per condition) than it returns.
int main(int argc, char **argv)
{
    return verify::run(argc, argv, {9});
}
