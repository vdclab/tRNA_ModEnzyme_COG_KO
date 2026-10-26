#!/usr/bin/env python3
"""
Figure: KO prevalence vs copy-number-when-present are uncorrelated.

Scatter of the 74 tRNA-modification KOs (x = prevalence, the fraction of
1,322 genomes with >=1 copy; y = mean copy number among genomes where the
KO is present). Points are colored by whether the KO ever reaches >=2
copies. The OLS fit and its 95% confidence band visualize the absence of
a relationship; r and rho are reported in the manuscript text, not on the
figure.

Input : /mnt/user-uploads/KO_count.txt
Output: fig3_dosage_vs_prevalence.svg/.png and
        fig3_dosage_vs_prevalence_points.tsv (plotted values)
        in /mnt/results/ko_copy_prevalence/

Usage: python3 script_fig_dosage_prevalence.py
Requires: pandas, numpy, matplotlib, scipy
"""

import os
import shutil

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as sps

KO_COUNT_FILE = "KO_count.txt"
OUTPUT_DIR = "ko_copy_prevalence"
WORK_DIR = "script_out"

BLUE, ORANGE = "#0279EE", "#FF9400"
LABEL_KOS = ["K00554", "K01011", "K06177", "K05812", "K04085"]  # extremes to label

plt.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"


def main():
    counts = pd.read_csv(KO_COUNT_FILE, sep="\t", index_col=0)
    pres = counts > 0
    df = pd.DataFrame({
        "prevalence": pres.mean(),
        "dosage_if_present": counts.where(pres).mean(),
        "max_copies": counts.max(),
    })
    df["ever_multicopy"] = df["max_copies"] >= 2

    # statistics (for the caption / text, not drawn on the figure)
    r, p_r = sps.pearsonr(df["prevalence"], df["dosage_if_present"])
    rho, p_rho = sps.spearmanr(df["prevalence"], df["dosage_if_present"])
    print(f"Pearson r = {r:+.3f} (p = {p_r:.3f})")
    print(f"Spearman rho = {rho:+.3f} (p = {p_rho:.3f})")

    # OLS fit + 95% confidence band
    slope, intercept, _, _, _ = sps.linregress(df["prevalence"], df["dosage_if_present"])
    xs = np.linspace(0, 1, 200)
    yhat = intercept + slope * xs
    resid = df["dosage_if_present"] - (intercept + slope * df["prevalence"])
    s_err = np.sqrt((resid**2).sum() / (len(df) - 2))
    xbar = df["prevalence"].mean()
    ssx = ((df["prevalence"] - xbar)**2).sum()
    ci = sps.t.ppf(0.975, len(df) - 2) * s_err * np.sqrt(1 / len(df) + (xs - xbar)**2 / ssx)

    # figure
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    for flag, color, label in [(False, BLUE, "single-copy KOs"),
                               (True, ORANGE, "KOs reaching \u22652 copies")]:
        sub = df[df["ever_multicopy"] == flag]
        ax.scatter(sub["prevalence"], sub["dosage_if_present"], s=38,
                   c=color, edgecolor="white", linewidth=0.6,
                   label=label, zorder=3)
    ax.plot(xs, yhat, color="grey", lw=1.4, zorder=2)
    ax.fill_between(xs, yhat - ci, yhat + ci, color="grey", alpha=0.18,
                    lw=0, zorder=1, label="OLS fit \u00b1 95% CI")

    for ko in LABEL_KOS:
        row = df.loc[ko]
        ax.annotate(ko, (row["prevalence"], row["dosage_if_present"]),
                    fontsize=7.5, xytext=(5, 4), textcoords="offset points")

    ax.set_xlabel("Prevalence (fraction of genomes with \u22651 copy)")
    ax.set_ylabel("Mean copy number | present")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(0.95, 1.60)
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    fig.tight_layout()

    # save
    os.makedirs(WORK_DIR, exist_ok=True)
    fig.savefig(f"{WORK_DIR}/fig3_dosage_vs_prevalence.svg")
    fig.savefig(f"{WORK_DIR}/fig3_dosage_vs_prevalence.png", dpi=200)
    df.round(4).to_csv(f"{WORK_DIR}/fig3_dosage_vs_prevalence_points.tsv", sep="\t")
    plt.close(fig)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for fname in ["fig3_dosage_vs_prevalence.svg", "fig3_dosage_vs_prevalence.png",
                  "fig3_dosage_vs_prevalence_points.tsv"]:
        dst = os.path.join(OUTPUT_DIR, fname)
        if os.path.exists(dst):
            os.remove(dst)
        shutil.copy(os.path.join(WORK_DIR, fname), dst)
    print(f"done — figure and point table in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
