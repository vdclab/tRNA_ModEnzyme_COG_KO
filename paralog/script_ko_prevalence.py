#!/usr/bin/env python3
"""
KO copy-number vs prevalence analysis across bacterial genomes.

Inputs
------
KO_count.txt       : tab-delimited matrix of KO copy counts (rows = genomes,
                     first column = genome ID, header = KO keys)
genome_phylum.txt  : tab-delimited genome -> phylum mapping
                     (header: genome, Phylum)

Outputs (written to OUTPUT_DIR)
------------------------------
ko_prevalence_copynumber.tsv   per-KO stats: prevalence, copy-number stats,
                               dosage conditional on presence, unweighted
                               phylum-mean prevalence, phylum heterogeneity
ko_prevalence_by_phylum.tsv    prevalence per KO x phylum (phyla with
                               >= MIN_PHYLUM_GENOMES genomes)
ko_sampling_distortion.tsv     per-KO sampling-weight distortion
                               (overall prevalence vs unweighted phylum mean)
fig1_overall_vs_phylum_mean.*  scatter of overall vs unweighted phylum-mean
                               prevalence, colored by phylum range
fig2_top_heterogeneous.*       heatmap of the most phylum-heterogeneous KOs

Usage
-----
python3 script_ko_prevalence.py

Requires: pandas, numpy, matplotlib
"""

import os
import shutil

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
INPUT_DIR = "."
KO_COUNT_FILE = os.path.join(INPUT_DIR, "KO_count.txt")
PHYLUM_FILE = os.path.join(INPUT_DIR, "genome_phylum.txt")

OUTPUT_DIR = "ko_copy_prevalence"
WORK_DIR = "script_out"

MIN_PHYLUM_GENOMES = 10   # phyla below this n are excluded from phylum tables
HET_TOP_N = 25            # KOs shown in the heterogeneity heatmap
HET_RANGE_THRESHOLD = 0.3 # phylum-range above which per-phylum discussion is advised

# Phylo palette
BLUE, ORANGE, GREEN, BEIGE = "#0279EE", "#FF9400", "#75A025", "#ECE9E2"

plt.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"


# --------------------------------------------------------------------------
# 1. Load inputs
# --------------------------------------------------------------------------
def load_inputs():
    counts = pd.read_csv(KO_COUNT_FILE, sep="\t", index_col=0)
    phylum = pd.read_csv(PHYLUM_FILE, sep="\t")

    # keep only genomes present in both files
    common = counts.index.intersection(phylum["genome"])
    dropped = counts.index.difference(phylum["genome"])
    if len(dropped):
        print(f"  note: {len(dropped)} genomes in KO_count.txt have no phylum "
              f"entry and were dropped")
    counts = counts.loc[common]
    pm = phylum.set_index("genome").loc[common, "Phylum"]

    print(f"  genomes: {counts.shape[0]}  KOs: {counts.shape[1]}  "
          f"phyla: {pm.nunique()}")
    return counts, pm


# --------------------------------------------------------------------------
# 2. Per-KO prevalence and copy-number statistics
# --------------------------------------------------------------------------
def per_ko_stats(counts, pm, big_phyla):
    pres = counts > 0

    # prevalence per KO x phylum (well-sampled phyla only)
    prev_phylum = pres.groupby(pm).mean().loc[big_phyla].T

    stats = pd.DataFrame({
        "n_present": pres.sum(),
        "prevalence": pres.mean(),
        "mean_copies_all": counts.mean(),
        "mean_copies_if_present": counts.where(pres).mean(),
        "max_copies": counts.max(),
        "min_copies_if_present": counts.where(pres).min(),
    })
    # unweighted mean of phylum prevalences (sampling-robust reference)
    stats["phylum_mean_prev_unweighted"] = prev_phylum.mean(axis=1)
    stats["phylum_prev_min"] = prev_phylum.min(axis=1)
    stats["phylum_prev_max"] = prev_phylum.max(axis=1)
    stats["phylum_prev_range"] = prev_phylum.max(axis=1) - prev_phylum.min(axis=1)
    stats["phylum_prev_sd"] = prev_phylum.std(axis=1)
    # sampling-weight distortion: overall (sample-weighted) minus unweighted
    stats["sampling_distortion"] = stats["prevalence"] - stats["phylum_mean_prev_unweighted"]
    # practical flag: per-phylum discussion advised?
    stats["phylum_specific"] = stats["phylum_prev_range"] > HET_RANGE_THRESHOLD

    # sanity check: mean copies == prevalence * dosage|present
    recon = stats["prevalence"] * stats["mean_copies_if_present"]
    assert np.allclose(recon, stats["mean_copies_all"], atol=1e-8), \
        "decomposition check failed"

    return stats.sort_values("prevalence", ascending=False), prev_phylum


# --------------------------------------------------------------------------
# 3. Figures
# --------------------------------------------------------------------------
def fig1_scatter(stats):
    """Overall (sample-weighted) vs unweighted phylum-mean prevalence."""
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    sc = ax.scatter(stats["phylum_mean_prev_unweighted"], stats["prevalence"],
                    c=stats["phylum_prev_range"], cmap="viridis",
                    s=40, edgecolor="white", linewidth=0.5, zorder=3)
    ax.plot([0, 1], [0, 1], ls="--", lw=1, color="grey", zorder=2)
    ax.set_xlabel("Unweighted mean of phylum prevalences")
    ax.set_ylabel("Overall prevalence (sample-weighted)")
    cb = fig.colorbar(sc, ax=ax, shrink=0.85)
    cb.set_label("Phylum prevalence range")
    # label the most distorted KOs
    dist = stats["sampling_distortion"].abs().nlargest(5).index
    for ko in dist:
        r = stats.loc[ko]
        ax.annotate(ko, (r["phylum_mean_prev_unweighted"], r["prevalence"]),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    fig.tight_layout()
    return fig


def fig2_heatmap(prev_phylum, stats):
    """Prevalence of the most phylum-heterogeneous KOs across phyla."""
    top = stats["phylum_prev_range"].nlargest(HET_TOP_N).index
    mat = prev_phylum.loc[top]
    order = mat.mean(axis=0).sort_values(ascending=False).index
    mat = mat[order]

    cmap = LinearSegmentedColormap.from_list("phylo", [BEIGE, BLUE])
    fig, ax = plt.subplots(figsize=(0.32 * mat.shape[1] + 2.2,
                                    0.28 * mat.shape[0] + 1.6))
    im = ax.imshow(mat.values, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=8)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("Prevalence")
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# 4. Main
# --------------------------------------------------------------------------
def main():
    os.makedirs(WORK_DIR, exist_ok=True)

    print("[1/4] Loading inputs")
    counts, pm = load_inputs()

    phylum_sizes = pm.value_counts()
    big_phyla = phylum_sizes[phylum_sizes >= MIN_PHYLUM_GENOMES].index.tolist()
    print(f"[2/4] {len(big_phyla)} phyla with >= {MIN_PHYLUM_GENOMES} genomes "
          f"({phylum_sizes[phylum_sizes >= MIN_PHYLUM_GENOMES].sum()} of "
          f"{len(pm)} genomes)")

    print("[3/4] Computing per-KO statistics")
    stats, prev_phylum = per_ko_stats(counts, pm, big_phyla)

    n_phylum_specific = stats["phylum_specific"].sum()
    print(f"  {n_phylum_specific}/{len(stats)} KOs have phylum prevalence "
          f"range > {HET_RANGE_THRESHOLD} (per-phylum discussion advised)")
    print(f"  max sampling distortion: {stats['sampling_distortion'].abs().max():.3f} "
          f"({stats['sampling_distortion'].abs().idxmax()})")

    print("[4/4] Writing outputs")
    stats.round(4).to_csv(f"{WORK_DIR}/ko_prevalence_copynumber.tsv", sep="\t")
    prev_phylum.round(4).to_csv(f"{WORK_DIR}/ko_prevalence_by_phylum.tsv", sep="\t")
    stats[["prevalence", "phylum_mean_prev_unweighted",
           "sampling_distortion"]].round(4).to_csv(
        f"{WORK_DIR}/ko_sampling_distortion.tsv", sep="\t")

    for name, fig in [("fig1_overall_vs_phylum_mean", fig1_scatter(stats)),
                      ("fig2_top_heterogeneous",
                       fig2_heatmap(prev_phylum, stats))]:
        fig.savefig(f"{WORK_DIR}/{name}.svg")
        fig.savefig(f"{WORK_DIR}/{name}.png", dpi=200)
        plt.close(fig)

    # sync to results (S3-backed mount forbids overwriting: remove first)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for fname in os.listdir(WORK_DIR):
        dst = os.path.join(OUTPUT_DIR, fname)
        if os.path.exists(dst):
            os.remove(dst)
        shutil.copy(os.path.join(WORK_DIR, fname), dst)

    print(f"done — all outputs in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
