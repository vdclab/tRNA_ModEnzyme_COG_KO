#!/usr/bin/env python3
"""
OLS regression equation for fig3 (dosage-when-present ~ prevalence).

Computes slope and intercept two ways:
  1. manually via the least-squares normal equations:
       slope     = sum((x-xbar)(y-ybar)) / sum((x-xbar)^2)
       intercept = ybar - slope * xbar
  2. via scipy.stats.linregress
and prints the fitted equation y = a + b*x with standard errors,
R^2, and its p-value.

Input : /mnt/user-uploads/KO_count.txt
Output: printed equation + ols_equation_dosage_prevalence.tsv
        in /mnt/results/ko_copy_prevalence/

Usage: python3 script_ols_equation.py
"""

import os
import shutil

import numpy as np
import pandas as pd
from scipy import stats as sps

KO_COUNT_FILE = "KO_count.txt"
OUTPUT_DIR = "ko_copy_prevalence"
WORK_DIR = "script_out"


def main():
    counts = pd.read_csv(KO_COUNT_FILE, sep="\t", index_col=0)
    pres = counts > 0
    x = pres.mean().values
    y = counts.where(pres).mean().values
    n = len(x)

    # 1) manual normal equations
    xbar, ybar = x.mean(), y.mean()
    sxy = ((x - xbar) * (y - ybar)).sum()
    sxx = ((x - xbar) ** 2).sum()
    slope = sxy / sxx
    intercept = ybar - slope * xbar

    # residual standard error and SEs of coefficients
    resid = y - (intercept + slope * x)
    s_err = np.sqrt((resid ** 2).sum() / (n - 2))
    se_slope = s_err / np.sqrt(sxx)
    se_intercept = s_err * np.sqrt(1 / n + xbar ** 2 / sxx)

    # 2) scipy reference
    lr = sps.linregress(x, y)
    assert np.isclose(lr.slope, slope) and np.isclose(lr.intercept, intercept)

    r, p_r = sps.pearsonr(x, y)

    print(f"n = {n}")
    print(f"manual  : slope = {slope:+.6f} (SE {se_slope:.6f}), "
          f"intercept = {intercept:.6f} (SE {se_intercept:.6f})")
    print(f"scipy   : slope = {lr.slope:+.6f} (SE {lr.stderr:.6f}), "
          f"intercept = {lr.intercept:.6f} (SE {lr.intercept_stderr:.6f})")
    print(f"\nOLS equation:  dosage_if_present = {intercept:.4f} "
          f"{slope:+.4f} x prevalence")
    print(f"  R^2 = {r**2:.5f}   p(slope) = {lr.pvalue:.3f}")

    out = pd.DataFrame([{
        "n": n,
        "intercept": intercept,
        "intercept_se": se_intercept,
        "slope": slope,
        "slope_se": se_slope,
        "r_squared": r ** 2,
        "p_value_slope": lr.pvalue,
        "equation": f"dosage_if_present = {intercept:.4f} {slope:+.4f} * prevalence",
    }])
    os.makedirs(WORK_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    out.round(6).to_csv(f"{WORK_DIR}/ols_equation_dosage_prevalence.tsv",
                        sep="\t", index=False)
    dst = os.path.join(OUTPUT_DIR, "ols_equation_dosage_prevalence.tsv")
    if os.path.exists(dst):
        os.remove(dst)
    shutil.copy(f"{WORK_DIR}/ols_equation_dosage_prevalence.tsv", dst)
    print(f"\nsaved: {dst}")


if __name__ == "__main__":
    main()
