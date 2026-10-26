#!/usr/bin/env python3
"""
Heatmaps of mean copy number per genome among genomes that carry the key,
COG (left panel) and KO (right panel), by phylum.

    cell = sum(copies over genomes in phylum with >= 1 copy)
           / (number of genomes in phylum with >= 1 copy)
    i.e. the phylum size is NOT the denominator; genomes with 0 copies are ignored.
    Cells with no carrier genome are left empty (light grey).
    Each cell shows its value to 1 decimal place.

    x = phyla (input-table order: largest first, "Other" last)
    y = keys, highest overall mean copy number among carriers at the top
    Cells are square. Both panels have the same height and share one colour scale.

Input is the by-phylum tables written by count_distribution.py.

Usage:
    python plot_mean_copy_heatmap.py COG_count_distribution_by_phylum.tsv \
                                     KO_count_distribution_by_phylum.tsv [outdir]

Outputs:
    COG_mean_copy_by_phylum.tsv
    KO_mean_copy_by_phylum.tsv
    COG_KO_mean_copy_heatmap.png
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FixedLocator, FixedFormatter

LOG_COLOUR = True       # log colour scale (most values are close to 1); False = linear
PANEL_HEIGHT = 14.0     # inches, height of each heatmap


def mean_copy_among_carriers(path):
    t = pd.read_csv(path, sep="\t", dtype={"Phylum": str})
    key = t.columns[0]
    count_cols = [c for c in t.columns if c not in (key, "Phylum", "N_genomes")]
    copies = np.array([int(c) for c in count_cols])
    phyla = list(dict.fromkeys(t["Phylum"]))
    sizes = t.groupby("Phylum")["N_genomes"].first()

    n = t[count_cols].to_numpy(float)                    # genomes with exactly k copies
    t["carriers"] = n.sum(axis=1)
    t["total_copies"] = (n * copies).sum(axis=1)

    carriers = t.pivot(index=key, columns="Phylum", values="carriers")[phyla]
    total = t.pivot(index=key, columns="Phylum", values="total_copies")[phyla]
    mean = total / carriers.where(carriers > 0)          # NaN where no carrier

    overall = total.sum(axis=1) / carriers.sum(axis=1).where(lambda s: s > 0)
    mean = mean.loc[overall.sort_values(ascending=False).index]
    mean.index.name = key
    return key, mean, sizes


def main(cog_file, ko_file, outdir="."):
    os.makedirs(outdir, exist_ok=True)
    panels = [mean_copy_among_carriers(cog_file), mean_copy_among_carriers(ko_file)]
    for key, mean, _ in panels:
        mean.round(4).to_csv(os.path.join(outdir, f"{key}_mean_copy_by_phylum.tsv"), sep="\t")

    # Shared colour range
    vmin = min(np.nanmin(m.to_numpy()) for _, m, _ in panels)
    vmax = max(np.nanmax(m.to_numpy()) for _, m, _ in panels)
    norm = LogNorm(vmin=max(1.0, vmin), vmax=vmax) if LOG_COLOUR else plt.Normalize(vmin, vmax)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e6e6e6")

    # Same height for both panels with square cells -> each panel's width is
    # height * (columns / rows); cell size differs slightly if row counts differ.
    widths = [PANEL_HEIGHT * m.shape[1] / m.shape[0] for _, m, _ in panels]
    label_space = 1.4                                    # inches for key labels on each side
    gap = 2.4                                            # inches between panels
    cbar_space = 1.6                                     # inches right of the KO labels
    fig_w = sum(widths) + 2 * label_space + gap + cbar_space + 0.3
    top, bottom = 1.0, 2.2                               # inches: title / phylum labels
    fig_h = PANEL_HEIGHT + top + bottom
    fig = plt.figure(figsize=(fig_w, fig_h))

    lefts = [label_space, label_space + widths[0] + gap]
    for i, ((key, mean, sizes), w, x0) in enumerate(zip(panels, widths, lefts)):
        ax = fig.add_axes((x0 / fig_w, bottom / fig_h, w / fig_w, PANEL_HEIGHT / fig_h))
        im = ax.imshow(np.ma.masked_invalid(mean.to_numpy()), cmap=cmap, norm=norm,
                       aspect="equal", interpolation="nearest")
        nr, nc = mean.shape
        # Cell values to 1 decimal place; font scaled to cell size, colour by background
        cell_pt = PANEL_HEIGHT / nr * 72
        fs = max(2.5, min(7.0, cell_pt * 0.36))
        vals = mean.to_numpy()
        for r in range(nr):
            for c in range(nc):
                v = vals[r, c]
                if np.isnan(v):
                    continue
                rgb = cmap(norm(v))[:3]
                lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
                ax.text(c, r, f"{v:.1f}", ha="center", va="center", fontsize=fs,
                        color="black" if lum > 0.5 else "white")
        ax.set_xticks(range(nc))
        ax.set_xticklabels([f"{p} ({sizes[p]})" for p in mean.columns],
                           rotation=90, fontsize=6.5)
        ax.set_yticks(range(nr))
        ax.set_yticklabels(mean.index, fontsize=max(4, min(7, 450 / nr)))
        if i == 1:                                       # KO labels on the outer (right) side
            ax.yaxis.tick_right()
        ax.set_xticks(np.arange(-0.5, nc, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, nr, 1), minor=True)
        ax.grid(which="minor", color="white", lw=0.3)
        ax.tick_params(which="minor", length=0)
        ax.set_title(f"{key} ({nr})", fontsize=11)

    # Vertical colour key on the right, centred on the panels
    cbar_h = PANEL_HEIGHT * 0.45
    cbar_x = lefts[1] + widths[1] + label_space + 0.2
    cax = fig.add_axes((cbar_x / fig_w, (bottom + (PANEL_HEIGHT - cbar_h) / 2) / fig_h,
                        0.22 / fig_w, cbar_h / fig_h))
    cb = fig.colorbar(im, cax=cax, orientation="vertical")
    if LOG_COLOUR:
        ticks = [t for t in [1, 1.5, 2, 3, 5, 10, 20, 50] if norm.vmin <= t <= norm.vmax]
        cb.ax.yaxis.set_major_locator(FixedLocator(ticks))
        cb.ax.yaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
        cb.ax.yaxis.set_minor_locator(FixedLocator([]))
    cb.ax.tick_params(labelsize=8)
    cb.set_label("Mean copies per genome\namong genomes with \u2265 1 copy"
                 + ("\n(log scale)" if LOG_COLOUR else "") + "; grey = no carrier", fontsize=9)
    fig.suptitle("Mean copy number of COGs and KOs by phylum (carrier genomes only)",
                 y=1 - 0.35 / fig_h, fontsize=13)

    out = os.path.join(outdir, "COG_KO_mean_copy_heatmap.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Colour range {norm.vmin:.2f}-{norm.vmax:.2f}; wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    main(*sys.argv[1:])
