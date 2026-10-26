#!/usr/bin/env python3
"""Step 11 (v2): Per-trait KO-KO within-pathway phi grids.

For each selected trait, render a grid of phi-coefficient heatmaps (one
panel per pathway with >= 3 KOs) from that trait's complete-case subset.
A shared color key is placed in the bottom-right corner of the grid.
Fonts slightly larger than v1. Saves .png + .svg in each trait's folder.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm

mpl.rcParams.update({
    "font.family": ["Liberation Sans", "Arimo", "DejaVu Sans"],
    "svg.fonttype": "none",
    "figure.dpi": 150,
})
BASE = Path("/mnt/results/ko_trait_analysis_v2")

for trait in ["Gram_Stain", "Motility", "Oxygen_Req.", "Growth_Temp_Category"]:
    csv_path = BASE / trait / f"ko_ko_within_pathway_{trait}.csv"
    if not csv_path.exists():
        print("missing, skip:", csv_path)
        continue
    df = pd.read_csv(csv_path)

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
    nrow = int(np.ceil(max(len(pathways), 1) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.6 * nrow))
    axes = np.atleast_2d(axes)
    flat = axes.ravel()
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    cmap = mpl.colormaps["RdBu_r"]

    n_plotted = 0
    for ax, pathway in zip(flat, pathways):
        m, q = mats[pathway]
        annot = pd.DataFrame("", index=m.index, columns=m.columns)
        annot[q == 0] = "**"
        annot[q == 1] = "*"
        sns.heatmap(m, cmap="RdBu_r", center=0, vmin=-1, vmax=1, ax=ax,
                    annot=annot, fmt="", annot_kws={"fontsize": 10},
                    cbar=False, square=True)
        ax.set_title(f"{pathway} (n={len(m)})", fontsize=13)
        ax.tick_params(labelsize=9)
        n_plotted += 1

    # shared color key in the bottom-right corner: use the last empty panel
    # if available, otherwise a small inset axes
    empty = [ax for ax in flat[n_plotted:]]
    if empty:
        cax = empty[-1]
        cax.axis("off")
        cax_pos = cax.get_position()
        cbar_ax = fig.add_axes([cax_pos.x0 + 0.02, cax_pos.y0 + 0.12,
                                cax_pos.width - 0.06, 0.035])
    else:
        cbar_ax = fig.add_axes([0.78, 0.03, 0.16, 0.02])
    cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm,
                                   orientation="horizontal")
    cb.set_label("phi coefficient", fontsize=11)
    cb.ax.tick_params(labelsize=10)

    fig.suptitle(f"KO-KO phi within pathway, {trait} subset "
                 "(* q<0.05, ** q<0.01, phylogenetic)", fontsize=14)
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(BASE / trait / f"ko_ko_phi_grid_{trait}.svg", bbox_inches="tight")
    fig.savefig(BASE / trait / f"ko_ko_phi_grid_{trait}.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote", BASE / trait / f"ko_ko_phi_grid_{trait}.png")
print("Step 11 done.")
