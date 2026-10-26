#!/usr/bin/env python3
"""Step 9: Figures for pairwise KO-KO within-pathway correlations.
Grid of phi-coefficient heatmaps, one panel per pathway with >= 3 KOs,
cells annotated with * (phylo q<0.05) / ** (q<0.01).
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

mpl.rcParams.update({
    "font.family": ["Liberation Sans", "Arimo", "DejaVu Sans"],
    "svg.fonttype": "none",
    "figure.dpi": 150,
})
BASE = Path("/mnt/results/ko_trait_analysis_v2/02_within_pathway")
df = pd.read_csv(BASE / "ko_ko_within_pathway_stats.csv")

# build per-pathway phi matrices
mats = {}
for pathway, g in df.groupby("pathway"):
    kos = sorted(set(g["KO_1"]) | set(g["KO_2"]))
    if len(kos) < 3:
        continue
    m = pd.DataFrame(np.eye(len(kos)), index=kos, columns=kos)
    q = pd.DataFrame(1.0, index=kos, columns=kos)
    for _, r in g.iterrows():
        m.loc[r.KO_1, r.KO_2] = m.loc[r.KO_2, r.KO_1] = r.phi
        q.loc[r.KO_1, r.KO_2] = q.loc[r.KO_2, r.KO_1] = (
            0 if r.q_phylo < 0.01 else 1 if r.q_phylo < 0.05 else 2)
    mats[pathway] = (m, q)

pathways = sorted(mats, key=lambda p: -len(mats[p][0]))
ncol = 3
nrow = int(np.ceil(len(pathways) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.9 * nrow))
axes = np.atleast_2d(axes)
for ax, pathway in zip(axes.ravel(), pathways):
    m, q = mats[pathway]
    annot = pd.DataFrame("", index=m.index, columns=m.columns)
    annot[q == 0] = "**"
    annot[q == 1] = "*"
    sns.heatmap(m, cmap="RdBu_r", center=0, vmin=-1, vmax=1, ax=ax,
                annot=annot, fmt="", annot_kws={"fontsize": 7},
                cbar=False, square=True)
    ax.set_title(f"{pathway} (n={len(m)})", fontsize=9)
    ax.tick_params(labelsize=6)
for ax in axes.ravel()[len(pathways):]:
    ax.axis("off")
fig.suptitle("KO-KO phi within pathway (* q<0.05, ** q<0.01, phylogenetic)",
             fontsize=11)
fig.tight_layout()
fig.savefig(BASE / "ko_ko_phi_grid.svg", bbox_inches="tight")
fig.savefig(BASE / "ko_ko_phi_grid.png", bbox_inches="tight")
print("wrote", BASE / "ko_ko_phi_grid.png")

# summary bar: fraction of pairs significant per pathway
summ = df.groupby("pathway").agg(
    n_pairs=("phi", "size"),
    n_sig_phylo=("q_phylo", lambda s: int((s < 0.05).sum())),
    mean_phi=("phi", "mean")).sort_values("mean_phi", ascending=False)
summ.to_csv(BASE / "ko_ko_pathway_summary.csv")
fig, ax = plt.subplots(figsize=(7, 4))
ax.barh(np.arange(len(summ)), summ["n_sig_phylo"] / summ["n_pairs"],
        color="#0279EE", height=0.7)
ax.set_yticks(np.arange(len(summ)))
ax.set_yticklabels([f"{p} ({s}/{n})" for p, s, n in
                    zip(summ.index, summ.n_sig_phylo, summ.n_pairs)], fontsize=8)
ax.set_xlabel("fraction of KO pairs significant (phylo q<0.05)")
ax.set_title("Within-pathway KO co-occurrence per pathway", fontsize=10)
fig.tight_layout()
fig.savefig(BASE / "ko_ko_fraction_sig.svg", bbox_inches="tight")
fig.savefig(BASE / "ko_ko_fraction_sig.png", bbox_inches="tight")
print("wrote", BASE / "ko_ko_fraction_sig.png")
