#!/usr/bin/env python3
"""Step 13: Per-trait pathway-pathway correlation heatmaps.

For each selected trait, render the 33x33 phylogeny-corrected correlation
matrix (r_phylo) as a clustered heatmap with a color key, using the SAME
pathway order across traits (from the all-genome analysis) for
comparability. Cells annotated * (q<0.05) / ** (q<0.01). Saves .png + .svg
in each trait's folder next to the source CSVs.
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

# fixed pathway order from the all-genome analysis (comparable across traits)
ref = pd.read_csv(BASE / "01_pathway_pathway" / "r_phylo_matrix.csv", index_col=0)
order = list(ref.columns[leaves_list(linkage(ref.values, method="average"))])

for trait in ["Gram_Stain", "Motility", "Oxygen_Req.", "Growth_Temp_Category"]:
    r_path = BASE / trait / f"pathway_pathway_r_matrix_{trait}.csv"
    q_path = BASE / trait / f"pathway_pathway_q_matrix_{trait}.csv"
    if not r_path.exists():
        print("missing, skip:", r_path)
        continue
    R = pd.read_csv(r_path, index_col=0).loc[order, order]
    Q = pd.read_csv(q_path, index_col=0).loc[order, order]

    annot = pd.DataFrame("", index=order, columns=order)
    annot[Q < 0.01] = "**"
    annot[(Q >= 0.01) & (Q < 0.05)] = "*"

    fig, ax = plt.subplots(figsize=(13, 10.5))
    sns.heatmap(R, cmap="RdBu_r", center=0, vmin=-1, vmax=1, ax=ax,
                cbar_kws={"label": "phylogeny-corrected r", "shrink": 0.7},
                annot=annot, fmt="", annot_kws={"fontsize": 8},
                xticklabels=True, yticklabels=True,
                linewidths=0.3, linecolor="0.95")
    ax.tick_params(labelsize=9)
    ax.set_title(f"Pathway x pathway correlation, {trait} subset "
                 "(* q<0.05, ** q<0.01, phylogenetic)", fontsize=13)
    fig.tight_layout()
    fig.savefig(BASE / trait / f"pathway_pathway_heatmap_{trait}.svg",
                bbox_inches="tight")
    fig.savefig(BASE / trait / f"pathway_pathway_heatmap_{trait}.png",
                bbox_inches="tight")
    plt.close(fig)
    print("wrote", BASE / trait / f"pathway_pathway_heatmap_{trait}.png")
print("Step 13 done.")
