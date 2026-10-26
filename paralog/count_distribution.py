#!/usr/bin/env python3
"""
Distribution of per-genome counts for each key (COG or KO), overall and by phylum.

Usage:
    python count_distribution.py COUNT_MATRIX GENOME_PHYLUM KEY_LABEL OUTDIR [MIN_GENOMES]

    COUNT_MATRIX   tab-delimited; first column genome IDs, header = keys
    GENOME_PHYLUM  tab-delimited; header "genome<TAB>Phylum"
    KEY_LABEL      label for the key column, e.g. COG or KO
    OUTDIR         output directory
    MIN_GENOMES    phyla with fewer genomes are pooled as "Other" (default 10).
                   Phylum names are used exactly as written: names that differ only
                   by a suffix (e.g. Bacillota, Bacillota_A) are never merged.

Outputs (prefix = KEY_LABEL):
    <KEY>_count_distribution.tsv            KEY, 1, 2, ..., max  (number of genomes)
    <KEY>_count_distribution_by_phylum.tsv  KEY, Phylum, N_genomes, 1, 2, ..., max
    <KEY>_phylum_groups.tsv                 original phylum -> plotted group
    <KEY>_count_distribution_heatmap.png
    <KEY>_count_distribution_by_phylum_heatmap.png
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def read_counts(path):
    m = pd.read_csv(path, sep="\t", index_col=0, dtype=str)
    m.columns = m.columns.str.strip()
    m.index = m.index.astype(str).str.strip()
    m = m.apply(pd.to_numeric, errors="raise").fillna(0).astype(int)
    if m.index.duplicated().any():
        sys.exit(f"Error: duplicated genome IDs in {path}")
    return m


def read_phylum(path):
    p = pd.read_csv(path, sep="\t", dtype=str)
    p.columns = p.columns.str.strip()
    p = p.apply(lambda s: s.str.strip())
    return p.set_index(p.columns[0])[p.columns[1]]


def distribution(m, max_count):
    """keys x (1..max_count): number of genomes with exactly that count."""
    cols = range(1, max_count + 1)
    return pd.DataFrame({k: (m == k).sum(axis=0) for k in cols}, index=m.columns)


def draw_heatmap(ax, data, norm, cmap, show_ylabels=True, title=None, ylab_size=5):
    arr = np.ma.masked_less_equal(data.to_numpy(dtype=float), 0)   # 0 -> white
    im = ax.imshow(arr, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    ax.set_xticks(range(data.shape[1]), data.columns, fontsize=6)
    if show_ylabels:
        ax.set_yticks(range(data.shape[0]), data.index, fontsize=ylab_size)
    else:
        ax.set_yticks([])
    ax.set_facecolor("white")
    if title:
        ax.set_title(title, fontsize=8)
    return im


def main(count_file, phylum_file, key, outdir, min_genomes=10):
    os.makedirs(outdir, exist_ok=True)
    m = read_counts(count_file)
    phylum = read_phylum(phylum_file)

    missing = m.index.difference(phylum.index)
    if len(missing):
        print(f"Warning: {len(missing)} genome(s) without phylum -> 'Unknown'", file=sys.stderr)
    ph = phylum.reindex(m.index).fillna("Unknown")

    max_count = int(m.to_numpy().max())
    print(f"{key}: {m.shape[0]} genomes x {m.shape[1]} keys, max count = {max_count}",
          file=sys.stderr)

    # ---- 1) overall distribution -----------------------------------------
    dist = distribution(m, max_count)
    dist.index.name = key
    dist.to_csv(os.path.join(outdir, f"{key}_count_distribution.tsv"), sep="\t")

    # ---- 2) by phylum ------------------------------------------------------
    sizes = ph.value_counts()
    group = ph.where(ph.map(sizes) >= min_genomes, "Other")
    gsizes = group.value_counts()
    order = [g for g in gsizes.index if g != "Other"] + (["Other"] if "Other" in gsizes else [])

    pd.DataFrame({"Phylum": sizes.index, "N_genomes": sizes.values,
                  "Group": [g if n >= min_genomes else "Other" for g, n in sizes.items()]}) \
        .to_csv(os.path.join(outdir, f"{key}_phylum_groups.tsv"), sep="\t", index=False)
    print(f"{len(order)} phylum groups (min {min_genomes} genomes; "
          f"{(sizes < min_genomes).sum()} small phyla pooled as 'Other' "
          f"with {gsizes.get('Other', 0)} genomes)", file=sys.stderr)

    blocks = {}
    for g in order:
        d = distribution(m[group == g], max_count)
        d.insert(0, "N_genomes", int(gsizes[g]))
        d.insert(0, "Phylum", g)
        blocks[g] = d
    by_ph = pd.concat(blocks.values())
    by_ph.index.name = key
    by_ph.to_csv(os.path.join(outdir, f"{key}_count_distribution_by_phylum.tsv"), sep="\t")

    # ---- 3) heatmaps -------------------------------------------------------
    # Rows ordered by mean count across all genomes (highest first), same in all panels
    row_order = m.mean(axis=0).sort_values(ascending=False).index
    cmap = plt.get_cmap("viridis").copy()
    nkeys = len(row_order)
    row_h = max(0.09, min(0.25, 10 / nkeys))        # inches per row, kept readable

    # 3a. overall: number of genomes
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * max_count + 2.5), nkeys * row_h + 1.5))
    im = draw_heatmap(ax, dist.loc[row_order], LogNorm(1, max(1, dist.to_numpy().max())), cmap)
    ax.set(xlabel=f"Copies of {key} per genome", ylabel=key,
           title=f"Distribution of {key} copy numbers across {m.shape[0]} genomes\n"
                 f"(white = no genome; genomes with 0 copies not shown)")
    fig.colorbar(im, ax=ax, label="Number of genomes (log scale)", shrink=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"{key}_count_distribution_heatmap.png"), dpi=300)
    plt.close(fig)

    # 3b. by phylum: % of genomes in that phylum, shared colour scale
    ncol = 5
    nrow = int(np.ceil(len(order) / ncol))
    panel_w = max(2.6, 0.13 * max_count + 1.0)
    facet_row_h = min(row_h, 4.5 / nkeys)           # keep each panel <= ~4.5 in tall
    fig, axes = plt.subplots(nrow, ncol, figsize=(panel_w * ncol + 1.5, (nkeys * facet_row_h + 0.9) * nrow),
                             squeeze=False)
    norm = LogNorm(0.1, 100)
    for i, g in enumerate(order):
        ax = axes[i // ncol][i % ncol]
        pct = blocks[g].drop(columns=["Phylum", "N_genomes"]).loc[row_order] / gsizes[g] * 100
        im = draw_heatmap(ax, pct, norm, cmap, show_ylabels=(i % ncol == 0),
                          title=f"{g} (n={gsizes[g]})", ylab_size=3.5)
    for j in range(len(order), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.supxlabel(f"Copies of {key} per genome", fontsize=10)
    fig.supylabel(key, fontsize=10)
    fig.suptitle(f"Distribution of {key} copy numbers by phylum "
                 f"(% of genomes in each phylum; white = none)", fontsize=12)
    fig.tight_layout(rect=(0.01, 0.01, 0.93, 0.98))
    cax = fig.add_axes((0.945, 0.35, 0.012, 0.3))
    fig.colorbar(im, cax=cax, label="% of genomes in phylum (log scale)")
    fig.savefig(os.path.join(outdir, f"{key}_count_distribution_by_phylum_heatmap.png"), dpi=300)
    plt.close(fig)

    print(f"{key}: outputs written to {outdir}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    a = sys.argv[1:5]
    main(*a, min_genomes=int(sys.argv[5]) if len(sys.argv) == 6 else 10)
