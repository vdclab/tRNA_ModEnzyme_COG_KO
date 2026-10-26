#!/usr/bin/env python3
"""
Side-by-side heatmaps of the multi-copy percentage by phylum.

    Left panel  = COG, right panel = KO
    x = phylum groups (largest first, "Other" last)
    y = keys, highest overall multi-copy % at the top
    cell = 100 x (genomes with >= 2 copies) / (genomes with >= 1 copy), within the phylum
    Grey cells = key absent from every genome of that phylum (percentage undefined).
    Cells are square, both panels have the same height and share one colour scale (0-100%).

Input are the by-phylum tables written by count_distribution.py.

Usage:
    python plot_multicopy_heatmap.py COG_by_phylum.tsv KO_by_phylum.tsv [outdir] [--min-present N]

    --min-present N  grey out cells where fewer than N genomes carry the key
                     (default 1 = only cells with no carrier are grey)

Outputs:
    COG_multicopy_pct_by_phylum.tsv, KO_multicopy_pct_by_phylum.tsv
    COG_KO_multicopy_heatmap.png
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def multicopy_table(path, min_present):
    t = pd.read_csv(path, sep="\t", dtype={"Phylum": str})
    key = t.columns[0]
    count_cols = [c for c in t.columns if c not in (key, "Phylum", "N_genomes")]
    multi_cols = [c for c in count_cols if int(c) >= 2]
    phyla = list(dict.fromkeys(t["Phylum"]))
    sizes = t.groupby("Phylum")["N_genomes"].first()

    t["present"] = t[count_cols].sum(axis=1)               # genomes with >= 1 copy
    t["multi"] = t[multi_cols].sum(axis=1)                 # genomes with >= 2 copies
    present = t.pivot(index=key, columns="Phylum", values="present")[phyla]
    multi = t.pivot(index=key, columns="Phylum", values="multi")[phyla]
    pct = (multi / present * 100).where(present >= max(min_present, 1))

    # Row order: multi-copy % over all genomes, highest first
    overall = multi.sum(axis=1) / present.sum(axis=1).replace(0, np.nan) * 100
    order = overall.sort_values(ascending=False, na_position="last").index
    pct, present = pct.loc[order], present.loc[order]
    pct.columns.name = present.columns.name = None
    return key, pct, present, sizes


def draw(ax, m, key, sizes, norm, cmap):
    arr = np.ma.masked_invalid(m.to_numpy(dtype=float))
    im = ax.imshow(arr, cmap=cmap, norm=norm, aspect="equal", interpolation="nearest")
    ax.set_xticks(range(m.shape[1]),
                  [f"{p} ({sizes[p]})" for p in m.columns], rotation=90, fontsize=7)
    ax.set_yticks(range(m.shape[0]), m.index, fontsize=6)
    ax.set_xticks(np.arange(-0.5, m.shape[1]), minor=True)
    ax.set_yticks(np.arange(-0.5, m.shape[0]), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4)
    ax.tick_params(which="minor", length=0)
    ax.set_title(f"{key} (n={m.shape[0]})", fontsize=11, pad=8)
    return im


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cog_table")
    ap.add_argument("ko_table")
    ap.add_argument("outdir", nargs="?", default=".")
    ap.add_argument("--min-present", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    panels = [multicopy_table(a.cog_table, a.min_present),
              multicopy_table(a.ko_table, a.min_present)]
    for key, pct, present, _ in panels:
        pct.round(2).rename_axis(key).to_csv(
            os.path.join(a.outdir, f"{key}_multicopy_pct_by_phylum.tsv"), sep="\t", na_rep="NA")

    norm = Normalize(0, 100)                               # shared colour range
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("lightgrey")

    # ---- explicit layout in inches: equal panel heights, square cells -------
    H = 14.0
    left, gap, cbar_gap, cbar_w, right = 1.0, 1.1, 0.35, 0.18, 0.9
    top, bottom = 0.7, 2.0
    widths = [H * pct.shape[1] / pct.shape[0] for _, pct, _, _ in panels]
    W = left + widths[0] + gap + widths[1] + cbar_gap + cbar_w + right
    FH = bottom + H + top
    fig = plt.figure(figsize=(W, FH))

    x = left
    for (key, pct, _, sizes), w in zip(panels, widths):
        ax = fig.add_axes((x / W, bottom / FH, w / W, H / FH))
        im = draw(ax, pct, key, sizes, norm, cmap)
        x += w + gap
    x = x - gap + cbar_gap
    cax = fig.add_axes((x / W, (bottom + H * 0.3) / FH, cbar_w / W, H * 0.4 / FH))
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("% of carrier genomes with ≥2 copies")

    grey = ("absent from phylum" if a.min_present <= 1
            else f"fewer than {a.min_present} carrier genomes")
    fig.suptitle(f"Multi-copy genomes among genomes carrying each key, by phylum "
                 f"(grey = {grey})", y=1 - 0.2 / FH, fontsize=12)
    out = os.path.join(a.outdir, "COG_KO_multicopy_heatmap.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
