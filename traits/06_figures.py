#!/usr/bin/env python3
"""Step 6: Figures for the KO-pathway vs trait analysis.

Reads the stats CSVs produced by the R steps and renders:
  - per-trait: forest plot of phylogenetic effects + significance bar heatmap
  - pathway x pathway: clustered heatmap of phylogeny-corrected correlations
  - within-pathway: forest plot of KO ~ LOO completeness
  - per-phylum: dot plot of subgroup effects
  - overview: mean pathway completeness per phylum
All figures saved as .svg (editable text) + .png.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage

mpl.rcParams.update({
    "font.family": ["Liberation Sans", "Arimo", "DejaVu Sans"],
    "svg.fonttype": "none",
    "figure.dpi": 150,
})
BASE = Path("/mnt/results/ko_trait_analysis_v2")
PALETTE = sns.color_palette("colorblind")


def save(fig, path: Path):
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", path.with_suffix(".png"))


def forest(df, est_col, lo_col, hi_col, label_col, title, out, xlabel):
    d = df.dropna(subset=[est_col]).sort_values(est_col)
    if d.empty:
        print("skip (empty):", title)
        return
    fig, ax = plt.subplots(figsize=(6, 0.32 * len(d) + 1.6))
    y = np.arange(len(d))
    ax.hlines(y, d[lo_col], d[hi_col], color="0.6", lw=1.2, zorder=1)
    colors = [PALETTE[0] if q < 0.05 else "0.55"
              for q in d.get("q_phylo", pd.Series([1] * len(d)))]
    ax.scatter(d[est_col], y, c=colors, s=28, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(d[label_col], fontsize=8)
    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=10)
    save(fig, out)


# ---------------------------------------------------------------- per-trait
TRAIT_LABEL = {
    "Gram_Stain": "Gram-positive vs negative (log OR)",
    "Motility": "Motile vs nonmotile (log OR)",
}
for trait, xlabel in TRAIT_LABEL.items():
    df = pd.read_csv(BASE / trait / "trait_pathway_stats.csv")
    df["lo"] = df["estimate"] - 1.96 * np.nan
    # phyloglm SEs are not in the CSV; plot estimates with significance color
    d = df.dropna(subset=["estimate"]).sort_values("estimate")
    fig, ax = plt.subplots(figsize=(7, 0.28 * len(d) + 1.6))
    colors = ["#0279EE" if q < 0.05 else "0.6" for q in d["q_phylo"]]
    ax.barh(np.arange(len(d)), d["estimate"], color=colors, height=0.7)
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(d["pathway"], fontsize=7)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_title(f"{trait}: pathway completeness effect (blue = q<0.05)",
                 fontsize=10)
    save(fig, BASE / trait / f"{trait}_forest")

    fig, ax = plt.subplots(figsize=(4.2, 6))
    pv = df.set_index("pathway")[["p_phylo"]].fillna(1.0)
    nlpv = -np.log10(pv.clip(lower=1e-12))
    sns.heatmap(nlpv, cmap="viridis", ax=ax, cbar_kws={"label": "-log10(p)"},
                yticklabels=True, xticklabels=True)
    ax.set_title(f"{trait}: phylogenetic test p-values", fontsize=10)
    ax.tick_params(axis="y", labelsize=6)
    save(fig, BASE / trait / f"{trait}_pvalue_heatmap")

# ------------------------------------------- categorical traits: F / epsilon^2
for trait in ["Oxygen_Req.", "Growth_Temp_Category"]:
    df = pd.read_csv(BASE / trait / "trait_pathway_stats.csv")
    d = df.dropna(subset=["F_phylo"]).sort_values("F_phylo")
    fig, ax = plt.subplots(figsize=(7, 0.28 * len(d) + 1.6))
    colors = ["#0279EE" if q < 0.05 else "0.6" for q in d["q_phylo"]]
    ax.barh(np.arange(len(d)), d["F_phylo"], color=colors, height=0.7)
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(d["pathway"], fontsize=7)
    ax.set_xlabel("phylogenetic ANOVA F")
    ax.set_title(f"{trait}: pathway completeness across states (blue = q<0.05)",
                 fontsize=10)
    save(fig, BASE / trait / f"{trait}_forest")

# --------------------------------------------------- pathway x pathway heatmap
R = pd.read_csv(BASE / "01_pathway_pathway" / "r_phylo_matrix.csv", index_col=0)
Q = pd.read_csv(BASE / "01_pathway_pathway" / "q_phylo_matrix.csv", index_col=0)
order = np.r_[leaves_list(linkage(R.loc[:, R.columns].values, method="average")),
              []]
order = [c for c in R.columns[leaves_list(linkage(R.values, method="average"))]]
R = R.loc[order, order]
Q = Q.loc[order, order]
annot = pd.DataFrame("", index=R.index, columns=R.columns)
annot[Q < 0.05] = "*"
fig, ax = plt.subplots(figsize=(11, 9))
sns.heatmap(R, cmap="RdBu_r", center=0, vmin=-1, vmax=1, ax=ax,
            cbar_kws={"label": "phylogeny-corrected r"}, annot=annot, fmt="",
            annot_kws={"fontsize": 6},
            xticklabels=True, yticklabels=True)
ax.tick_params(labelsize=6)
ax.set_title("Pathway x pathway correlation (completeness, * = q<0.05)",
             fontsize=10)
save(fig, BASE / "01_pathway_pathway" / "pathway_pathway_heatmap")

# ------------------------------------------------------------- within-pathway
_wp = BASE / "02_within_pathway" / "within_pathway_stats.csv"
if not _wp.exists():
    wp = None
    print("within-pathway stats not found (skipped in v2); skipping that figure")
else:
    wp = pd.read_csv(_wp)
if wp is not None:
    wp["label"] = wp["KO"] + " (" + wp["pathway"] + ")"
    fig, ax = plt.subplots(figsize=(7, 4))
    d = wp.dropna(subset=["estimate"]).sort_values("estimate")
    colors = ["#0279EE" if q < 0.05 else "0.6" for q in d["q_phylo"]]
    ax.barh(np.arange(len(d)), d["estimate"], color=colors, height=0.7)
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(d["label"], fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("phylogenetic logistic coefficient (KO presence ~ LOO completeness)")
    ax.set_title("Within-pathway KO vs pathway (blue = q<0.05)", fontsize=10)
    save(fig, BASE / "02_within_pathway" / "within_pathway_forest")

# ---------------------------------------------------------------- per-phylum
pp = Path(BASE / "03_per_phylum" / "per_phylum_stats.csv")
if pp.exists():
    df = pd.read_csv(pp)
    if len(df):
        df["pair"] = df["trait"] + "\n" + df["pathway"]
        n_pair = df["pair"].nunique()
        fig, ax = plt.subplots(figsize=(0.5 * n_pair + 4, 0.4 * df["phylum"].nunique() + 3))
        pivot = df.pivot_table(index="phylum", columns="pair", values="estimate",
                               aggfunc="first")
        sns.heatmap(pivot, cmap="RdBu_r", center=0, ax=ax,
                    cbar_kws={"label": "effect estimate"},
                    linewidths=0.5, linecolor="0.9")
        ax.set_title("Per-phylum subgroup effects (significant global pairs)",
                     fontsize=10)
        save(fig, BASE / "03_per_phylum" / "per_phylum_effects")

# ------------------------------------------------- phylum completeness overview
comp = pd.read_csv(BASE / "00_processed" / "pathway_completeness.csv", index_col=0)
traits = pd.read_csv(BASE / "00_processed" / "traits_selected.csv", index_col=0)
counts = comp.join(traits["Phylum"]).groupby("Phylum").size()
mean_comp = comp.join(traits["Phylum"]).groupby("Phylum").mean()
mean_comp = mean_comp.loc[counts.sort_values(ascending=False).index]
fig, ax = plt.subplots(figsize=(11, 6))
sns.heatmap(mean_comp.T, cmap="mako", ax=ax,
            cbar_kws={"label": "mean completeness"})
ax.set_title("Mean pathway completeness per phylum (phyla ordered by n genomes)",
             fontsize=10)
ax.tick_params(axis="x", labelsize=6, rotation=90)
ax.tick_params(axis="y", labelsize=7)
save(fig, BASE / "00_processed" / "phylum_completeness_overview")

print("All figures done.")
