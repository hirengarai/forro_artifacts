"""
complexity_5round.py

Data and time complexity of the divide-and-conquer PNB phase of the
5-round attack on Forro (Section IV-B), in the multi-list
form of S. Dey, "Advancing the Idea of Probabilistic Neutral Bits: First Key
Recovery Attack on 7.5 Round ChaCha", IEEE Trans. Inf. Theory 70(8), 2024
(https://doi.org/10.1109/TIT.2024.3389874).

The significant bits are split into k groups, each guessed against its own list,
so the time is a sum of 2^{m_i} N terms rather than one 2^m N. The PNB counts
and the conditional correlations come from Table VI.

Gives N = 2^28.51 and C = 2^107.97.

by Hiren

Run: python3 complexity_5round.py    (settings under __main__)
"""

import math
from statistics import NormalDist

# ----------------------------- Utility -----------------------------

def bias_product(biases):
    """
    Return the product of a list of biases.
    Example: biases=[2**-4.5, 2**-1.2] -> 2**(-5.7)
    decimal value also works
    """
    prod = 1.0
    for b in biases:
        prod *= b
    return prod


def compute_stage_N(alpha, fwd_eps, bwd_eps, p_nd=None):
    """
    Compute the required data complexity N for a single stage.

        epsilon = fwd_eps * bwd_eps
        one_minus_eps_sq = 1 - epsilon^2
        numerator = sqrt(alpha * ln 4) - constant * sqrt(1 - epsilon^2)
        N = (numerator / epsilon)^2

    If p_nd is not provided, defaults to NormalDist(mu=0, sigma=1).inv_cdf(0.0013)
    """

    if p_nd is None:
        constant = NormalDist(mu=0, sigma=1).inv_cdf(0.0013)
    else:
        constant = NormalDist(mu=0, sigma=1).inv_cdf(p_nd)

    epsilon = fwd_eps * bwd_eps
    print(f"The product of fwd and bwd correlations:~2^{{{math.log2(epsilon):.2f}}}.")
    one_minus_eps_sq = 1.0 - (epsilon * epsilon)
    numerator = math.sqrt(alpha * math.log(4.0)) - constant * math.sqrt(one_minus_eps_sq)
    N = (numerator / epsilon) ** 2

    return N


def compute_C(m_list, N, dim_g, R, r, key_size, alpha):
    """
    Compute and print log2 of each term in:

        C = sum_i 2^{m_i} * N
          + 2^{dim(g_new)} * N * (k - 1) / (2^{11.17} * (R - r))
          + 2^{|K| - alpha}
          + 2^{|K| - dim(g_new)}
    """
    k = len(m_list)

    term1 = sum(2 ** m for m in m_list) * N
    term2 = (2 ** dim_g) * N * (k - 1) / (2 ** 11.17 * (R - r))
    term3 = 2 ** (key_size - alpha)
    term4 = 2 ** (key_size - dim_g)

    print(f"\nTerm1 (sum 2^m_i * N)  :~2^{{{math.log2(term1):.2f}}}.")
    print(f"Term2 (list build)     :~2^{{{math.log2(term2):.2f}}}.")
    print(f"Term3 (2^(|K|-alpha))  :~2^{{{key_size - alpha:.2f}}}.")
    print(f"Term4 (2^(|K|-dim g))  :~2^{{{key_size - dim_g:.2f}}}.")

    C = term1 + term2 + term3 + term4
    print(f"\nFinal time complexity:~2^{{{math.log2(C):.2f}}}.")
    if C >= 2 ** key_size:
        print(f"[!] exceeds 2^{{{key_size}}} brute force.")
    print()
    return C


# ----------------------------- Forro 5 rounds -----------------------------
if __name__ == "__main__":
    key_size = 196   # key bits left for this phase
    init_pnb_count = 107
    dim_g = key_size - init_pnb_count
    p_nd = 0.0013
    alpha = 95.46

    pnb_per_bit = [42, 55, 27, 50]
    bwd_corr_per_pnb = [0.93, 1, 0.88, 1]

    total_round, dist_round = 5, 3
    fwd_eps = pow(2, -7.92)
    bwd_corr = 0.22

    bwd_eps = bwd_corr * bias_product(bwd_corr_per_pnb)
    print(f"fwd corr: {fwd_eps:.5f} ~ 2^{{{math.log2(fwd_eps):.2f}}}, "
          f"bwd corr: {bwd_eps:.5f} ~ 2^{{{math.log2(bwd_eps):.2f}}}.")

    # Each PNB group must have exactly one associated backward correlation.
    assert len(pnb_per_bit) == len(bwd_corr_per_pnb), (
        f"Length mismatch: pnb_per_bit has {len(pnb_per_bit)} entries, "
        f"bwd_corr_per_pnb has {len(bwd_corr_per_pnb)}"
    )

    sig_per_bit = []
    for count in pnb_per_bit:
        sig_per_bit.append(dim_g - count)

    N = compute_stage_N(alpha, fwd_eps, bwd_eps, p_nd)
    print(f"The init. data complexity:~2^{{{math.log2(N):.2f}}}.")

    # Compute C using that N
    compute_C(sig_per_bit, N, dim_g, total_round, dist_round, key_size, alpha)
