#!/usr/bin/env python3
"""
Explore COG-KO concordance across genomes and phyla.

Usage:
    python cog_ko_explore.py COG.txt KO.txt key.txt genome_phylum.txt outdir [min_genomes]

min_genomes (default 10): phyla with fewer genomes are pooled as "Other" in plots.

Outputs (in outdir):
  genome_summary.tsv   per-genome counts, agreement, Jaccard, within-phylum outlier flag
  pair_summary.tsv     per-pair counts, agreement, Jaccard, phi correlation
  pair_by_phylum_jaccard.tsv
  pair_by_phylum_jaccard_sorted.tsv   same, rows in fig4 heatmap order
  fig1_genome_jaccard_by_phylum.png
  fig2_genome_ranked.png
  fig3_pair_distributions.png
  fig4_pair_phylum_heatmap.png
  fig5_pcoa_genomes.png
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def read_matrix(path):
    m = pd.read_csv(path, sep="\t", index_col=0, dtype=str)
    m.columns = m.columns.str.strip()
    m.index = m.index.astype(str).str.strip()
    return m.apply(pd.to_numeric).astype(np.int8)


def safe_div(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.divide(a, b, out=np.full_like(a, np.nan), where=b > 0)


def main(cog_file, ko_file, key_file, phylum_file, outdir, min_genomes=10):
    os.makedirs(outdir, exist_ok=True)
    cog, ko = read_matrix(cog_file), read_matrix(ko_file)
    genomes = cog.index[cog.index.isin(ko.index)]
    cog, ko = cog.loc[genomes], ko.loc[genomes]

    keys = pd.read_csv(key_file, sep="\t", dtype=str)
    keys.columns = keys.columns.str.strip()
    keys = keys[["COG", "KO"]].apply(lambda s: s.str.strip()).dropna().drop_duplicates()
    keys = keys[keys["COG"].isin(cog.columns) & keys["KO"].isin(ko.columns)].reset_index(drop=True)
    pair_names = keys["COG"] + "|" + keys["KO"]

    phy = pd.read_csv(phylum_file, sep="\t", dtype=str)
    phy = phy.apply(lambda s: s.str.strip())
    phy = phy.set_index(phy.columns[0])[phy.columns[1]]
    phylum = phy.reindex(genomes).fillna("Unknown")
    n_missing = (phylum == "Unknown").sum()
    if n_missing:
        print(f"Warning: {n_missing} genome(s) have no phylum -> 'Unknown'", file=sys.stderr)

    # ---- 4-state comparison (genomes x pairs) -----------------------------
    C = cog[keys["COG"]].to_numpy(bool)
    K = ko[keys["KO"]].to_numpy(bool)
    both, none = C & K, ~C & ~K
    cog_only, ko_only = C & ~K, ~C & K
    agree = (C == K)

    # ---- Per-genome summary ----------------------------------------------
    g = pd.DataFrame({
        "phylum": phylum.values,
        "both_present": both.sum(1), "both_absent": none.sum(1),
        "COG_only": cog_only.sum(1), "KO_only": ko_only.sum(1),
    }, index=genomes)
    g["agreement"] = agree.mean(1)
    g["jaccard"] = safe_div(g["both_present"], g["both_present"] + g["COG_only"] + g["KO_only"])
    g["COG_only_share"] = safe_div(g["COG_only"], g["COG_only"] + g["KO_only"])

    # Robust z-score of Jaccard within phylum (median / MAD); global for tiny phyla
    def robust_z(x):
        med = x.median()
        mad = (x - med).abs().median() * 1.4826
        return (x - med) / mad if mad > 0 else x * 0.0
    counts = g["phylum"].map(g["phylum"].value_counts())
    g["robust_z_within_phylum"] = np.where(
        counts >= min_genomes,
        g.groupby("phylum")["jaccard"].transform(robust_z),
        robust_z(g["jaccard"]))
    g["outlier"] = g["robust_z_within_phylum"].abs() > 3.5
    g.sort_values("jaccard").to_csv(os.path.join(outdir, "genome_summary.tsv"), sep="\t")

    # ---- Per-pair summary --------------------------------------------------
    n = len(genomes)
    a, b, c, d = both.sum(0), cog_only.sum(0), ko_only.sum(0), none.sum(0)
    phi = safe_div(a * d - b * c, np.sqrt((a + b) * (c + d) * (a + c) * (b + d)))
    p = pd.DataFrame({
        "COG": keys["COG"], "KO": keys["KO"],
        "COG_prevalence": C.mean(0), "KO_prevalence": K.mean(0),
        "both_present": a, "COG_only": b, "KO_only": c, "both_absent": d,
        "agreement": agree.mean(0), "jaccard": safe_div(a, a + b + c), "phi": phi,
    }, index=pair_names)
    p.sort_values("jaccard").to_csv(os.path.join(outdir, "pair_summary.tsv"), sep="\t")

    # ---- Plot grouping: pool small phyla ---------------------------------
    vc = g["phylum"].value_counts()
    g["phylum_plot"] = g["phylum"].where(g["phylum"].map(vc) >= min_genomes, "Other")
    order = (g.groupby("phylum_plot")["jaccard"].median().sort_values().index.tolist())
    labels = {k: f"{k} (n={v})" for k, v in g["phylum_plot"].value_counts().items()}
    palette = dict(zip(order, sns.color_palette("tab20", len(order))))

    # Fig 1: per-genome Jaccard by phylum
    fig, ax = plt.subplots(figsize=(8, 0.4 * len(order) + 2))
    sns.boxplot(data=g, y="phylum_plot", x="jaccard", order=order, color="white",
                fliersize=0, ax=ax)
    sns.stripplot(data=g[~g["outlier"]], y="phylum_plot", x="jaccard", order=order,
                  size=2.5, alpha=0.5, color="grey", ax=ax)
    sns.stripplot(data=g[g["outlier"]], y="phylum_plot", x="jaccard", order=order,
                  size=4, color="red", ax=ax)
    ax.set_yticks(range(len(order)), [labels[o] for o in order])
    ax.set(xlabel="COG-KO Jaccard (both present / either present)", ylabel="",
           title="Per-genome concordance by phylum\n(red = within-phylum outlier)")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig1_genome_jaccard_by_phylum.png"), dpi=200)
    plt.close(fig)

    # Fig 2: genomes ranked by Jaccard, coloured by phylum
    gs = g.sort_values("jaccard").reset_index()
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for ph in order:
        s = gs[gs["phylum_plot"] == ph]
        ax.scatter(s.index, s["jaccard"], s=6, color=palette[ph], label=labels[ph])
    out = gs[gs["outlier"]]
    ax.scatter(out.index, out["jaccard"], s=30, facecolors="none", edgecolors="red", lw=1)
    for _, r in out.head(15).iterrows():
        ax.annotate(r[gs.columns[0]], (r.name, r["jaccard"]), fontsize=6,
                    xytext=(3, -3), textcoords="offset points")
    ax.set(xlabel="Genome rank", ylabel="COG-KO Jaccard",
           title="Genomes ranked by concordance (red circles = outliers)")
    ax.legend(fontsize=6, markerscale=2, ncol=2, loc="lower right", frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig2_genome_ranked.png"), dpi=200)
    plt.close(fig)

    # Fig 3: per-pair distributions
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    axes[0].hist(p["jaccard"].dropna(), bins=50, color="steelblue")
    axes[0].set(xlabel="Jaccard", ylabel="Number of COG|KO pairs", title="Per-pair Jaccard")
    axes[1].hist(p["phi"].dropna(), bins=50, color="darkorange")
    axes[1].set(xlabel="Phi correlation", title="Per-pair phi")
    sc = axes[2].scatter(p["COG_prevalence"], p["KO_prevalence"], s=4, c=p["jaccard"],
                         cmap="viridis", vmin=0, vmax=1)
    axes[2].plot([0, 1], [0, 1], "k--", lw=0.7)
    axes[2].set(xlabel="COG prevalence", ylabel="KO prevalence", title="Prevalence per pair")
    fig.colorbar(sc, ax=axes[2], label="Jaccard")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig3_pair_distributions.png"), dpi=200)
    plt.close(fig)

    # Fig 4: pair x phylum Jaccard heatmap for the most phylum-variable pairs
    grp = g["phylum_plot"].to_numpy()
    rows = {}
    for ph in order:
        m = grp == ph
        rows[ph] = safe_div(both[m].sum(0), (both | cog_only | ko_only)[m].sum(0))
    pp = pd.DataFrame(rows, index=pair_names)
    pp.to_csv(os.path.join(outdir, "pair_by_phylum_jaccard.tsv"), sep="\t")
    if len(pp) >= 2:
        # Cluster rows on Jaccard; NaN (pair absent from a whole phylum) -> -0.1 for
        # clustering only, and shown in grey via the mask.
        show_labels = len(pp) <= 150
        cg = sns.clustermap(pp.fillna(-0.1), cmap="RdYlBu", vmin=0, vmax=1,
                            col_cluster=False, figsize=(10, 12),       # height:width = 1.2:1
                            yticklabels=show_labels, xticklabels=True,
                            mask=pp.isna(), dendrogram_ratio=(0.12, 0.04),
                            cbar_pos=(0.45, 0.975, 0.3, 0.01),
                            cbar_kws={"label": "Jaccard", "orientation": "horizontal"})
        cg.ax_heatmap.set_facecolor("lightgrey")
        cg.ax_cbar.xaxis.set_ticks_position("top")
        cg.ax_cbar.xaxis.set_label_position("top")
        cg.ax_heatmap.tick_params(axis="x", labelsize=7)
        cg.ax_heatmap.tick_params(axis="y", labelsize=4)
        cg.ax_heatmap.set_xticklabels([labels[t.get_text()] for t in cg.ax_heatmap.get_xticklabels()],
                                      rotation=90)
        cg.ax_heatmap.set_xlabel("")
        cg.ax_heatmap.set_ylabel(f"COG|KO pairs (n={len(pp)})" +
                                 ("" if show_labels else ", labels hidden"))
        cg.fig.suptitle("COG-KO concordance (Jaccard) of all pairs across phyla\n"
                        "(grey = pair not found in that phylum)", y=1.065)
        # Same table, rows in heatmap order (top to bottom)
        row_order = cg.dendrogram_row.reordered_ind
        sorted_pp = pp.iloc[row_order].copy()
        sorted_pp.insert(0, "heatmap_row", np.arange(1, len(sorted_pp) + 1))
        sorted_pp.index.name = "COG|KO"
        sorted_pp.to_csv(os.path.join(outdir, "pair_by_phylum_jaccard_sorted.tsv"), sep="\t")
        cg.savefig(os.path.join(outdir, "fig4_pair_phylum_heatmap.png"), dpi=300,
                   bbox_inches="tight")
        plt.close(cg.fig)

    # Fig 5: PCoA of genomes on discordance profile (0 agree, -1 KO only, +1 COG only)
    X = cog_only.astype(np.int8) - ko_only.astype(np.int8)
    X = X[:, X.any(0)].astype(float)
    if X.shape[1] >= 2:
        sq = (X ** 2).sum(1)
        D2 = sq[:, None] + sq[None, :] - 2 * X @ X.T          # squared Euclidean
        J = np.eye(n) - 1.0 / n
        B = -0.5 * J @ D2 @ J
        vals, vecs = np.linalg.eigh(B)
        idx = np.argsort(vals)[::-1][:2]
        coords = vecs[:, idx] * np.sqrt(np.clip(vals[idx], 0, None))
        expl = vals[idx] / vals[vals > 0].sum() * 100
        fig, ax = plt.subplots(figsize=(8, 6))
        for ph in order:
            m = grp == ph
            ax.scatter(coords[m, 0], coords[m, 1], s=8, alpha=0.7, color=palette[ph], label=labels[ph])
        ax.set(xlabel=f"PCo1 ({expl[0]:.1f}%)", ylabel=f"PCo2 ({expl[1]:.1f}%)",
               title="PCoA of genomes by COG-KO discordance pattern")
        ax.legend(fontsize=6, markerscale=2, bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
        fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig5_pcoa_genomes.png"), dpi=200)
        plt.close(fig)

    print(f"{n} genomes, {len(keys)} pairs, {g['outlier'].sum()} outlier genome(s). "
          f"Results in {outdir}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) not in (6, 7):
        sys.exit(__doc__)
    args = sys.argv[1:6]
    mg = int(sys.argv[6]) if len(sys.argv) == 7 else 10
    main(*args, min_genomes=mg)
