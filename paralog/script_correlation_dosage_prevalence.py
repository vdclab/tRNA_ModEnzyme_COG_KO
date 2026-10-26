#!/usr/bin/env python3
"""
Correlation between KO prevalence and copy-number-when-present across KOs:
explicit calculation of Pearson r, its p-value, and Spearman rho.

The script computes each statistic three ways:
  1. manually from the defining formulas (transparency),
  2. via scipy.stats (the standard implementation),
  3. via a permutation test (robust p-value, recommended at n = 74).

Input : /mnt/user-uploads/KO_count.txt (copy counts, genomes x KO)
Output: printed comparison + correlation_dosage_prevalence.tsv
        in /mnt/results/ko_copy_prevalence/

Usage: python3 script_correlation_dosage_prevalence.py
Requires: pandas, numpy, scipy
"""

import os
import shutil

import numpy as np
import pandas as pd
from scipy import stats as sps

KO_COUNT_FILE = "KO_count.txt"
OUTPUT_DIR = "ko_copy_prevalence"
WORK_DIR = "script_out"
N_PERM = 100_000
SEED = 42


# ----------------------------------------------------------------------
# Load data and build the two per-KO vectors
# ----------------------------------------------------------------------
def load_vectors():
    counts = pd.read_csv(KO_COUNT_FILE, sep="\t", index_col=0)
    pres = counts > 0
    x = pres.mean().values                       # prevalence per KO
    y = counts.where(pres).mean().values         # dosage | present per KO
    return x, y


# ----------------------------------------------------------------------
# 1. Pearson r, manually
#    r = sum((x-mx)(y-my)) / sqrt(sum((x-mx)^2) * sum((y-my)^2))
#    t = r * sqrt(n-2) / sqrt(1-r^2)   -> two-sided p from t(n-2)
# ----------------------------------------------------------------------
def pearson_manual(x, y):
    n = len(x)
    xm, ym = x - x.mean(), y - y.mean()
    r = (xm * ym).sum() / np.sqrt((xm**2).sum() * (ym**2).sum())
    t = r * np.sqrt(n - 2) / np.sqrt(1 - r**2)
    p = 2 * sps.t.sf(np.abs(t), df=n - 2)
    return r, t, p


# ----------------------------------------------------------------------
# 2. Spearman rho, manually
#    rank both vectors (average ranks for ties), then Pearson on ranks.
#    p-value: t-approximation, same formula as above on rho.
# ----------------------------------------------------------------------
def rankdata_average(v):
    """Average ranks with ties, equivalent to scipy.stats.rankdata."""
    order = np.argsort(v, kind="mergesort")
    ranks = np.empty(len(v))
    ranks[order] = np.arange(1, len(v) + 1, dtype=float)
    # replace tied blocks by their average rank
    sv = v[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    return ranks


def spearman_manual(x, y):
    rx, ry = rankdata_average(x), rankdata_average(y)
    rho, t, p = pearson_manual(rx, ry)
    return rho, t, p


# ----------------------------------------------------------------------
# 3. Permutation p-values (two-sided): shuffle y, recompute statistic
# ----------------------------------------------------------------------
def permutation_p(x, y, stat_fn, n_perm=N_PERM, seed=SEED):
    rng = np.random.default_rng(seed)
    obs = stat_fn(x, y)[0]
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = stat_fn(x, rng.permutation(y))[0]
    p = (np.sum(np.abs(null) >= np.abs(obs)) + 1) / (n_perm + 1)
    return obs, p


# ----------------------------------------------------------------------
def main():
    x, y = load_vectors()
    n = len(x)
    print(f"n KOs = {n}")
    print(f"prevalence        : mean {x.mean():.3f}  range [{x.min():.3f}, {x.max():.3f}]")
    print(f"dosage | present  : mean {y.mean():.3f}  range [{y.min():.3f}, {y.max():.3f}]\n")

    # --- Pearson -------------------------------------------------------
    r_m, t_m, p_m = pearson_manual(x, y)
    r_sci, p_sci = sps.pearsonr(x, y)
    print("Pearson r (prevalence vs dosage|present)")
    print(f"  manual   : r = {r_m:+.4f}   t({n}-2) = {t_m:+.4f}   p = {p_m:.4f}")
    print(f"  scipy    : r = {r_sci:+.4f}   p = {p_sci:.4f}")
    r_perm, p_perm = permutation_p(x, y, pearson_manual)
    print(f"  perm ({N_PERM:,}): r = {r_perm:+.4f}   p = {p_perm:.4f}\n")

    # --- Spearman ------------------------------------------------------
    rho_m, t_m, p_m = spearman_manual(x, y)
    rho_sci, p_sci = sps.spearmanr(x, y)
    print("Spearman rho (prevalence vs dosage|present)")
    print(f"  manual   : rho = {rho_m:+.4f}   t-approx p = {p_m:.4f}")
    print(f"  scipy    : rho = {rho_sci:+.4f}   p = {p_sci:.4f}")
    rho_perm, p_perm = permutation_p(x, y, spearman_manual)
    print(f"  perm ({N_PERM:,}): rho = {rho_perm:+.4f}   p = {p_perm:.4f}\n")

    # --- save ----------------------------------------------------------
    out = pd.DataFrame({
        "measure": ["pearson_r", "spearman_rho"],
        "value": [r_m, rho_m],
        "p_t_approx": [pearson_manual(x, y)[2], spearman_manual(x, y)[2]],
        "p_permutation": [permutation_p(x, y, pearson_manual)[1],
                          permutation_p(x, y, spearman_manual)[1]],
        "n": [n, n],
    })
    os.makedirs(WORK_DIR, exist_ok=True)
    dst = os.path.join(OUTPUT_DIR, "correlation_dosage_prevalence.tsv")
    out.round(6).to_csv(f"{WORK_DIR}/correlation_dosage_prevalence.tsv", sep="\t", index=False)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    if os.path.exists(dst):
        os.remove(dst)
    shutil.copy(f"{WORK_DIR}/correlation_dosage_prevalence.tsv", dst)
    print(f"saved: {dst}")


if __name__ == "__main__":
    main()
