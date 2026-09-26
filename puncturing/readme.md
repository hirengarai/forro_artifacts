# Experimental verification of the puncturing attack

Two experiments on 5-round Forro, both from Appendix C.1 of Florez-Gutierrez
and Todo, [*Improved Cryptanalysis of ChaCha: Beating PNBs with Bit
Puncturing*](https://eprint.iacr.org/2025/437) (EUROCRYPT 2025). `rho2`
measures `rho^2 = <f, g>`; `rho4` checks that a pair factorises as
`p = c * c'`.

## Steps

```
xz -dk data/f_1_5round.punc.xz                 # 1.2 MiB -> 88 MiB, once [decompress, unzipping]
g++ -std=c++20 -O3 -march=native -pthread -o rho2 rho2.cpp
g++ -std=c++20 -O3 -march=native -pthread -o rho4 rho4.cpp
./rho2 data/f_1_5round.punc
./rho4 data/f_1_5round.punc
```

Needs a C++20 compiler and under a GiB of memory. Each returns 0 if the check
passes, 1 if it fails. As shipped they run 1024 keys and take hours; for a
quick look set `KEYS = 8` and `SAMPLES_PER_KEY = 1 << 9` in the Control Panel
and rebuild, where 8 keys scatter by about 5%. `-v` adds the header, progress,
controls and verdict; a trailing number fixes the seed (`... 777 -v`) and
reproduces bit for bit on any machine and thread count.

## Expected output

Seed 777, so these reproduce exactly. Cost is linear in `KEYS` times
`SAMPLES_PER_KEY` and inverse in cores - 17 and 45 core-minutes per key -
so the quick look above is a few minutes on a laptop.

```
$ ./rho2 data/f_1_5round.punc 777 -v
   102 of 1024
   204 of 1024
   306 of 1024
   408 of 1024
   510 of 1024
   612 of 1024
   714 of 1024
   816 of 1024
   918 of 1024
  1020 of 1024
 ──────────────────────────────────────────────────────────────
 ρ² claimed        2^-0.87
 ρ² measured       2^-0.86
 ‖g‖²              2^-0.86           a self-check on g alone; must equal ρ²
 elapsed           2.9 h             10513 s, 1024 keys x 4096 samples on 96 threads
 agreement         ratio 1.007    +2.00 sigma    +-0.0019 over 1024 keys
 variance expected 2^-12.86          ‖g‖² / samples per key
 variance observed 2^-8.07           across 1024 keys, 27.6x predicted
 wrong key         2^-5.47           +-0.0023    +0.62 sigma from the floor
 wrong-key floor   2^-5.56           g averaged over every wrong guess
 ──────────────────────────────────────────────────────────────
 result            agrees
 ──────────────────────────────────────────────────────────────
```

Without `-v` only the lines from `ρ² claimed` to `elapsed` are printed.

The right-key variance running well above prediction is the key dependence
reported for ChaCha in Appendix C.1; Forro shows it more strongly.

The pair check:

```
$ ./rho4 data/f_1_5round.punc 777 -v
   102 of 1024
   204 of 1024
   306 of 1024
   408 of 1024
   510 of 1024
   612 of 1024
   714 of 1024
   816 of 1024
   918 of 1024
  1020 of 1024
 ──────────────────────────────────────────────────────────────
 ρ⁴ predicted      2^-1.69           c · c', measured per key and averaged
 ρ⁴ measured       2^-1.69           the pair correlation
 elapsed           3.9 h             13864 s, 1024 keys x 4096 samples on 96 threads
 agreement         ratio 1.000    -0.18 sigma    +-0.0004 over 1024 keys
 variance expected 2^-13.73          ρ⁴ / samples per key
 variance observed 2^-7.71           across 1024 keys, 65.1x predicted
 mean c            0.552561          per key, the single-text correlation
 mean p            0.308904          per key, the pair
 p - c · c'         -0.000055        +-0.000121    -0.45 sigma from 0
 (ρ²)²             2^-1.73           -2.5% from the prediction; E[c²] ≠ (E[c])², which is why the test is per key
 wrong key         2^-7.39           reported, not a pass/fail
 ──────────────────────────────────────────────────────────────
 result            factorises
 ──────────────────────────────────────────────────────────────
```

Without `-v` only the lines from `ρ⁴ predicted` to `elapsed` are printed.

`f` is the true target bit, read from inside the encryption; `g` is the
punctured map, computed only from the keystream and a guessed key. The second
predicts with `c * c'` per key, never the key-averaged `(rho^2)^2`, which is
1.5% away by Jensen; `-v` prints it alongside.

The map holds 2,822,720 terms for target `v14[27]` after subround 18, reads
36 keystream bits and guesses 60 key bits.

```
sha256  05ade5f1e2d3385b3b27a5162d95071e1d2cb1bad468ea90e43c47a9c44a2828
        2140466505685c686b58f14bd046078ee25323ea80ba806717564ff3c90db43f  (.xz)
```
