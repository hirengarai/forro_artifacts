# Forró artifacts

This is the code and the extra material behind our paper *Two-Phase Key
Recovery of Forró: Bit Puncturing Meets Probabilistic Neutral Bits*. There is
nothing to install. You only need a C++ compiler that understands C++20 (g++ 10
or clang 10 will do) and Python 3.

There are three experiments, one per folder, plus a folder with some helping figures and
tables.

## syncopation: the 6.5- and 6.75-round attacks

Here we check the conditional backward correlation of the two syncopation
attacks. Each program reads its CPNB set from the text file next to it, works
out the restricted bits and the syncopated segments with Algorithm 1 of Wang et
al. (CRYPTO 2023), and then measures |eps'_a| as the median over random keys.

    cd syncopation
    make
    ./run65
    ./run675

`run65` is the 6.5-round attack with 113 CPNBs and should print about 0.0425.
`run675` is the 6.75-round attack with 80 CPNBs and should print about 0.00861.
Both programs first test the cipher implementation and stop if that test
fails. The last digit can move a little from one machine to another.

## puncturing: the 5-round attack

Two checks of the punctured map of the 5-round attack, in the spirit of
Appendix C.1 of Flórez-Gutiérrez and Todo (EUROCRYPT 2025). `rho2` measures the
puncturing correlation ρ² of the map. `rho4` checks that the correlation of a
pair of keystreams is the product of the two single correlations.

    cd puncturing
    xz -dk data/f_1_5round.punc.xz
    g++ -std=c++20 -O3 -march=native -pthread -o rho2 rho2.cpp
    g++ -std=c++20 -O3 -march=native -pthread -o rho4 rho4.cpp
    ./rho2 data/f_1_5round.punc 777 -v
    ./rho4 data/f_1_5round.punc 777 -v

The paper claims ρ² = 2^-0.87 and the first program measures 2^-0.86. The
second one gives 2^-1.69 for both the prediction and the measurement. Be
warned that each run takes a few hours, since it goes through 1024 keys. The
777 is the seed, so you get exactly the output listed in
[puncturing/readme.md](puncturing/readme.md).

## complexity: the PNB phase of the 5-, 5.5- and 6-round attacks

Three small scripts that evaluate the complexity formulas of the PNB phase with
the parameters from the paper.

    cd complexity
    python3 complexity_5round.py
    python3 complexity_5p5round.py
    python3 complexity_6round.py

Each script prints the data N and the time C of its phase: 2^28.51 and
2^107.97 for 5 rounds, 2^40.22 and 2^112.38 for 5.5 rounds, and 2^43.82 and
2^124.25 for 6 rounds.

## paper-supplement

Some useful figures and tables: the
trail figures and the mask-extension figures of the 5-, 5.5- and 6-round
attacks, the key-bit pictures of the 6.5- and 6.75-round attacks, the
per-component PNB sets and the 5.5-round trail tables. Every PDF has a caption
that explains it.
