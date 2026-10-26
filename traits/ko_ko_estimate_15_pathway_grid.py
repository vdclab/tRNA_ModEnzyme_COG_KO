#!/usr/bin/env python3
"""
KO-KO co-occurrence grids within pathway, Gram_Stain subset.

Each cell shows the phylogenetically corrected association between two KEGG
Orthology (KO) pairs from phylogenetic logistic regression (phylolm::phyloglm,
MPLE): the `estimate` column = log-odds coefficient. Positive (red) =
co-occurrence, negative (blue) = mutual exclusion, 0 (white) = no association.
Scale is clipped at +/-3 log-odds because the raw estimates are heavily
right-skewed (up to ~21).

Significance stars come from the FDR-adjusted phylogenetic p-value (q_phylo):
* q < 0.05, ** q < 0.01.

The diagonal (a KO paired with itself) is undefined on the log-odds scale and
is masked grey.

All 15 pathways in the input file are plotted; pathways with only 2 KOs
(one pair) appear as 2x2 grids.
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ----------------------------- configuration ------------------------------
INPUT_CSV = "/mnt/user-uploads/ko_ko_within_pathway_Gram_Stain.csv"
OUT_PNG = "/mnt/results/ko_ko_estimate_grid_Gram_Stain_all15.png"
OUT_SVG = "/mnt/results/ko_ko_estimate_grid_Gram_Stain_all15.svg"
VMAX = 3            # clip log-odds color scale at +/- VMAX
KO_FONTSIZE = 12    # KO tick-label font size (1.5 x the original 8)
STAR_FONTSIZE = 10  # significance marker font size
MIN_KOS = 2         # include pathways with at least this many KOs
# ---------------------------------------------------------------------------

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

df = pd.read_csv(INPUT_CSV)

# pathways ordered by number of KOs (descending), then alphabetically
sizes = df.groupby("pathway").apply(
    lambda g: len(set(g["KO_1"]) | set(g["KO_2"])), include_groups=False
)
pathways = sorted(sizes[sizes >= MIN_KOS].index, key=lambda p: (-sizes[p], p))

cmap = matplotlib.colormaps["RdBu_r"].copy()
cmap.set_bad("#eeeeee")  # masked diagonal
norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)

ncols = 3
nrows = int(np.ceil(len(pathways) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(17, 5.2 * nrows))
axes = np.atleast_2d(axes)

for idx, pw in enumerate(pathways):
    r, c = divmod(idx, ncols)
    ax = axes[r][c]
    sub = df[df["pathway"] == pw]
    kos = sorted(set(sub["KO_1"]) | set(sub["KO_2"]))
    n = len(kos)
    pos = {k: i for i, k in enumerate(kos)}

    M = np.full((n, n), np.nan)          # diagonal stays NaN -> grey
    star = np.full((n, n), "", dtype=object)
    for _, row in sub.iterrows():
        i, j = pos[row["KO_1"]], pos[row["KO_2"]]
        M[i, j] = M[j, i] = row["estimate"]
        q = row["q_phylo"]
        s = "**" if q < 0.01 else ("*" if q < 0.05 else "")
        star[i, j] = star[j, i] = s

    im = ax.imshow(M, cmap=cmap, norm=norm)
    for i in range(n):
        for j in range(n):
            if star[i, j]:
                v = M[i, j]
                dark = np.isnan(v) or abs(norm(v) - 0.5) > 0.25
                ax.text(j, i, star[i, j], ha="center", va="center",
                        fontsize=STAR_FONTSIZE,
                        color="white" if dark else "black")

    ax.set_xticks(range(n), kos, rotation=90, fontsize=KO_FONTSIZE)
    ax.set_yticks(range(n), kos, fontsize=KO_FONTSIZE)
    ax.set_title(f"{pw} (n={n})", fontsize=13)
    ax.set_xticks(np.arange(-0.5, n), minor=True)
    ax.set_yticks(np.arange(-0.5, n), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="both", length=0)

# hide unused axes (15 pathways fill a 5x3 grid exactly, so usually none)
for k in range(len(pathways), nrows * ncols):
    r, c = divmod(k, ncols)
    axes[r][c].axis("off")

# shared horizontal colorbar below the grid
cax = fig.add_axes([0.36, 0.045, 0.28, 0.008])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("phylogeny-corrected log-odds (phyloglm estimate), clipped at \u00b13",
             fontsize=10)

fig.suptitle("KO-KO phyloglm estimate within pathway, Gram_Stain subset\n"
             "(* q<0.05, ** q<0.01, phylogenetic; grey diagonal = self-pair, not defined)",
             fontsize=15)
fig.tight_layout(rect=[0, 0.06, 1, 0.96])

fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
fig.savefig(OUT_SVG, bbox_inches="tight")
print(f"saved: {OUT_PNG}")
print(f"saved: {OUT_SVG}")
