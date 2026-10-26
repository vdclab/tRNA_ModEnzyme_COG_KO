#!/usr/bin/env python3
"""
Correlation between tRNA-modification KO presence and genome size.

Inputs
------
KO_AP.txt                : tab-delimited binary matrix (0/1) of KO presence
                           (rows = genomes, first column = genome ID)
genome_phylum_size.txt   : tab-delimited genome -> Phylum, size_Mbp
                           (header: genome, Phylum, Total_KO_AP_Count, size_Mbp)

Outputs (written to OUTPUT_DIR)
------------------------------
ko_size_correlation.tsv        per-KO correlations with genome size
                               (Pearson/point-biserial r, Spearman rho,
                               p-values, BH-FDR q-values, prevalence)
sum_ko_size_correlation.tsv    genome-level: sum of 74 KOs vs size
                               (correlations + OLS regression)
correlation_linear_vs_log_long.tsv   linear vs log(size) correlations,
                               pooled and per phylum (tidy: one row per
                               scope x model; Pearson r, Spearman rho,
                               p-values, OLS R2, slope, AIC)
correlation_linear_vs_log_wide.tsv   side-by-side comparison with dAIC
                               and the better-fitting model per scope
ko_size_correlation_by_phylum.tsv  per-phylum correlation of KO-count sum
                               vs size (phyla with >= MIN_PHYLUM_GENOMES)
fig1_sum_vs_size.*             regression plot: KO-count sum vs genome size
fig2_per_ko_correlations.*     per-KO Pearson r bar chart with FDR significance
fig3_linear_vs_log_examples.*  per-phylum scatter with linear and log(size)
                               OLS fits for all well-sampled phyla

Usage
-----
python3 script_ko_size_correlation.py

Requires: pandas, numpy, scipy, matplotlib
"""

import os
import shutil

import numpy as np
import pandas as pd
from scipy import stats as sps
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
INPUT_DIR = "."
KO_AP_FILE = os.path.join(INPUT_DIR, "KO_AP.txt")
META_FILE = os.path.join(INPUT_DIR, "genome_phylum_size.txt")

OUTPUT_DIR = "ko_size_correlation"
WORK_DIR = "script_out"

MIN_PHYLUM_GENOMES = 10
FDR_ALPHA = 0.05

# Phylo palette
BLUE, ORANGE, GREEN, BEIGE = "#0279EE", "#FF9400", "#75A025", "#ECE9E2"

plt.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"


# --------------------------------------------------------------------------
# 1. Load inputs
# --------------------------------------------------------------------------
def load_inputs():
    ap = pd.read_csv(KO_AP_FILE, sep="\t", index_col=0)
    meta = pd.read_csv(META_FILE, sep="\t").set_index("genome")

    common = ap.index.intersection(meta.index)
    if len(common) < len(ap):
        print(f"  note: dropped {len(ap) - len(common)} genomes without metadata")
    ap = ap.loc[common]
    meta = meta.loc[common]
    print(f"  genomes: {ap.shape[0]}  KOs: {ap.shape[1]}  "
          f"phyla: {meta['Phylum'].nunique()}")
    return ap, meta


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR q-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


# --------------------------------------------------------------------------
# 2. Per-KO correlations with genome size
# --------------------------------------------------------------------------
def per_ko_correlations(ap, meta):
    size = meta["size_Mbp"]
    rows = []
    for ko in ap.columns:
        x = ap[ko].values
        pr = sps.pearsonr(x, size)    # point-biserial: binary vs continuous
        sr = sps.spearmanr(x, size)
        rows.append({
            "KO": ko,
            "prevalence": x.mean(),
            "pearson_r": pr[0], "pearson_p": pr[1],
            "spearman_rho": sr[0], "spearman_p": sr[1],
        })
    df = pd.DataFrame(rows).set_index("KO")
    df["pearson_q"] = bh_fdr(df["pearson_p"])
    df["spearman_q"] = bh_fdr(df["spearman_p"])
    df["significant_fdr"] = (df["pearson_q"] < FDR_ALPHA) | (df["spearman_q"] < FDR_ALPHA)
    return df.sort_values("pearson_r")


# --------------------------------------------------------------------------
# 3. Genome-level: sum of KOs vs size
# --------------------------------------------------------------------------
def sum_correlation(ap, meta):
    ko_sum = ap.sum(axis=1)
    size = meta["size_Mbp"]

    pr = sps.pearsonr(ko_sum, size)
    sr = sps.spearmanr(ko_sum, size)
    reg = sps.linregress(size, ko_sum)

    summary = pd.DataFrame([{
        "n_genomes": len(ko_sum),
        "n_kos": ap.shape[1],
        "pearson_r": pr[0], "pearson_p": pr[1],
        "spearman_rho": sr[0], "spearman_p": sr[1],
        "ols_slope": reg.slope, "ols_intercept": reg.intercept,
        "ols_r2": reg.rvalue ** 2,
        "ols_slope_p": reg.pvalue,
        "ols_stderr": reg.stderr,
        "mean_ko_count": ko_sum.mean(), "sd_ko_count": ko_sum.std(),
        "mean_size_Mbp": size.mean(), "sd_size_Mbp": size.std(),
    }])
    return ko_sum, summary


# --------------------------------------------------------------------------
# 4. Sampling-evenness check: per-phylum correlations (linear AND log)
# --------------------------------------------------------------------------
def per_phylum_correlations(ap, meta, ko_sum):
    """Per-phylum correlation of KO-count sum vs genome size, computed
    against raw size (linear) and log(size). Columns are suffixed
    _linear / _logsize; dAIC and better_model compare the two OLS fits."""
    size = meta["size_Mbp"]
    ph = meta["Phylum"]
    sizes = ph.value_counts()
    big = sizes[sizes >= MIN_PHYLUM_GENOMES].index

    rows = []
    for p in big:
        m = ph == p
        y = ko_sum[m]
        s = size[m]
        ls = np.log(s)
        pr_l, pr_g = sps.pearsonr(y, s), sps.pearsonr(y, ls)
        sr_l, sr_g = sps.spearmanr(y, s), sps.spearmanr(y, ls)
        reg_l, reg_g = sps.linregress(s, y), sps.linregress(ls, y)
        aic_l, aic_g = _gaussian_aic(s, y), _gaussian_aic(ls, y)
        rows.append({
            "Phylum": p, "n_genomes": int(m.sum()),
            "mean_size_Mbp": s.mean(), "sd_size_Mbp": s.std(),
            "mean_ko_count": y.mean(), "sd_ko_count": y.std(),
            # linear: count ~ size
            "pearson_r_linear": pr_l[0], "pearson_p_linear": pr_l[1],
            "spearman_rho_linear": sr_l[0], "spearman_p_linear": sr_l[1],
            "ols_r2_linear": reg_l.rvalue ** 2,
            "slope_linear": reg_l.slope, "aic_linear": aic_l,
            # log: count ~ log(size)
            "pearson_r_logsize": pr_g[0], "pearson_p_logsize": pr_g[1],
            "spearman_rho_logsize": sr_g[0], "spearman_p_logsize": sr_g[1],
            "ols_r2_logsize": reg_g.rvalue ** 2,
            "slope_logsize": reg_g.slope, "aic_logsize": aic_g,
            # comparison
            "dAIC_linear_minus_log": aic_l - aic_g,
            "better_model": "log_size" if aic_l > aic_g else "linear",
        })
    by_phylum = pd.DataFrame(rows).set_index("Phylum")

    # phylum-level confounding: do phylum means co-vary?
    pr_conf = sps.pearsonr(by_phylum["mean_size_Mbp"], by_phylum["mean_ko_count"])
    conf = pd.DataFrame([{
        "comparison": "phylum_means_size_vs_kocount",
        "n_phyla": len(by_phylum),
        "pearson_r": pr_conf[0], "pearson_p": pr_conf[1],
    }])
    return by_phylum, conf, big


# --------------------------------------------------------------------------
# 4b. Linear vs log(size) model comparison (pooled + per phylum)
# --------------------------------------------------------------------------
def _gaussian_aic(x, y):
    """AIC of a two-parameter OLS fit with Gaussian errors."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    X = np.column_stack([np.ones(len(x)), x])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    rss = float(((y - X @ b) ** 2).sum())
    n = len(y)
    ll = -n / 2 * (np.log(2 * np.pi) + np.log(rss / n) + 1)
    return 2 * 2 - 2 * ll  # k = 2 (slope + intercept); sigma is a function of rss


def correlation_comparison(ko_sum, meta, big_phyla):
    """Pearson/Spearman/AIC for count ~ size and count ~ log(size),
    pooled and per well-sampled phylum."""
    size = meta["size_Mbp"]
    ph = meta["Phylum"]

    scopes = [("POOLED", np.ones(len(ko_sum), dtype=bool))] + \
             [(p, (ph == p).values) for p in big_phyla]

    rows = []
    for scope, m in scopes:
        y = ko_sum[m].values
        s = size[m].values
        for model, x in [("linear", s), ("log_size", np.log(s))]:
            pr = sps.pearsonr(x, y)
            sr = sps.spearmanr(x, y)
            reg = sps.linregress(x, y)
            rows.append({
                "scope": scope, "n": int(m.sum()), "model": model,
                "pearson_r": pr[0], "pearson_p": pr[1],
                "spearman_rho": sr[0], "spearman_p": sr[1],
                "ols_r2": reg.rvalue ** 2,
                "ols_slope": reg.slope, "ols_intercept": reg.intercept,
                "aic": _gaussian_aic(x, y),
            })
    long = pd.DataFrame(rows)

    # wide side-by-side comparison
    wide = long.pivot(index="scope", columns="model",
                      values=["pearson_r", "pearson_p", "spearman_rho",
                              "spearman_p", "ols_r2", "aic"])
    wide.columns = [f"{stat}_{model}" for stat, model in wide.columns]
    wide["n"] = long.groupby("scope")["n"].first()
    wide["dAIC_linear_minus_log"] = wide["aic_linear"] - wide["aic_log_size"]
    wide["better_model"] = np.where(wide["dAIC_linear_minus_log"] > 0,
                                    "log_size", "linear")
    wide = wide.reset_index()
    return long, wide


def fig3_linear_vs_log(ko_sum, meta, big_phyla):
    """Grid of per-phylum scatters with linear and log(size) OLS fits."""
    size = meta["size_Mbp"]
    ph = meta["Phylum"]
    n = len(big_phyla)
    ncol = 5
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.9 * nrow))
    axes = np.atleast_2d(axes)

    for idx, p in enumerate(big_phyla):
        ax = axes[idx // ncol, idx % ncol]
        m = (ph == p).values
        s, y = size[m].values, ko_sum[m].values
        ax.scatter(s, y, s=10, alpha=0.5, color=BLUE, linewidths=0)
        xs = np.linspace(s.min(), s.max(), 100)
        reg_l = sps.linregress(s, y)
        reg_g = sps.linregress(np.log(s), y)
        ax.plot(xs, reg_l.intercept + reg_l.slope * xs, color=ORANGE,
                lw=1.6, label="linear")
        ax.plot(xs, reg_g.intercept + reg_g.slope * np.log(xs), color=GREEN,
                lw=1.6, ls="--", label="log(size)")
        ax.set_title(f"{p} (n={int(m.sum())})", fontsize=8)
        ax.tick_params(labelsize=7)
        if idx == 0:
            ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    for idx in range(n, nrow * ncol):
        axes[idx // ncol, idx % ncol].axis("off")
    for ax in axes[-1, :min(n, ncol)]:
        ax.set_xlabel("Genome size (Mbp)", fontsize=8)
    for r in range(nrow):
        axes[r, 0].set_ylabel("KO count", fontsize=8)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# 5. Figures
# --------------------------------------------------------------------------
def fig1_regression(ko_sum, meta, summary):
    size = meta["size_Mbp"]
    ph = meta["Phylum"]
    sizes = ph.value_counts()
    big = sizes[sizes >= MIN_PHYLUM_GENOMES].index
    col = ph.where(ph.isin(big), "other (<10 genomes)")

    palette = dict(zip(big, [BLUE, ORANGE, GREEN, "#B03A2E", "#7D3C98",
                             "#148F77", "#B7950B", "#5D6D7E", "#CA6F1E",
                             "#2E4053", "#943126", "#1F618D", "#884EA0",
                             "#A04000", "#117A65", "#797D7F", "#5499C7",
                             "#D68910"][:len(big)]))
    palette["other (<10 genomes)"] = "#D5D8DC"

    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    for p in list(big) + ["other (<10 genomes)"]:
        m = (col == p).values
        ax.scatter(size[m], ko_sum[m], s=14, alpha=0.65, linewidths=0,
                   color=palette[p], label=p, zorder=3)

    s = summary.iloc[0]
    xs = np.linspace(size.min(), size.max(), 100)
    ax.plot(xs, s["ols_intercept"] + s["ols_slope"] * xs, color="black",
            lw=1.8, zorder=4,
            label=(f"OLS: y = {s['ols_intercept']:.1f} + "
                   f"{s['ols_slope']:.2f}x  (R$^2$ = {s['ols_r2']:.2f})"))
    ax.set_xlabel("Genome size (Mbp)")
    ax.set_ylabel(f"Number of KOs present (of {summary.iloc[0]['n_kos']})")
    ax.legend(fontsize=6.5, frameon=False, loc="lower right", ncol=2)
    fig.tight_layout()
    return fig


def fig2_per_ko(df):
    d = df.sort_values("pearson_r")
    colors = [BLUE if r >= 0 else ORANGE for r in d["pearson_r"]]
    fig, ax = plt.subplots(figsize=(5.6, 9.5))
    ax.barh(d.index, d["pearson_r"], color=colors, height=0.75)
    for i, (ko, r) in enumerate(d["pearson_r"].items()):
        if d.loc[ko, "significant_fdr"]:
            star = "**" if min(d.loc[ko, "pearson_q"], d.loc[ko, "spearman_q"]) < 0.01 else "*"
            ax.text(r + (0.008 if r >= 0 else -0.008), i, star,
                    va="center", ha="left" if r >= 0 else "right", fontsize=8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Pearson r (KO presence vs genome size)")
    ax.set_xlim(-0.75, 0.75)
    ax.text(0.98, 0.02, f"* FDR < {FDR_ALPHA}   ** FDR < 0.01",
            transform=ax.transAxes, ha="right", fontsize=8, color="0.3")
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# 6. Main
# --------------------------------------------------------------------------
def main():
    os.makedirs(WORK_DIR, exist_ok=True)

    print("[1/5] Loading inputs")
    ap, meta = load_inputs()

    print("[2/5] Per-KO correlations with genome size")
    ko_corr = per_ko_correlations(ap, meta)
    nsig = ko_corr["significant_fdr"].sum()
    print(f"  {nsig}/{len(ko_corr)} KOs significant at FDR < {FDR_ALPHA}")
    print(f"  r range: [{ko_corr['pearson_r'].min():.3f}, "
          f"{ko_corr['pearson_r'].max():.3f}]")

    print("[3/5] Genome-level KO-count sum vs size")
    ko_sum, summary = sum_correlation(ap, meta)
    s = summary.iloc[0]
    print(f"  Pearson r = {s['pearson_r']:.3f} (p = {s['pearson_p']:.2e}), "
          f"Spearman rho = {s['spearman_rho']:.3f}, "
          f"OLS R2 = {s['ols_r2']:.3f}")

    print("[4/5] Per-phylum correlations, linear vs log(size)")
    by_phylum, conf, big_phyla = per_phylum_correlations(ap, meta, ko_sum)
    print(f"  {len(by_phylum)} phyla with >= {MIN_PHYLUM_GENOMES} genomes")
    print(f"  within-phylum Pearson r (linear): median "
          f"{by_phylum['pearson_r_linear'].median():.3f} "
          f"(range {by_phylum['pearson_r_linear'].min():.3f} to "
          f"{by_phylum['pearson_r_linear'].max():.3f})")
    print(f"  within-phylum Pearson r (log size): median "
          f"{by_phylum['pearson_r_logsize'].median():.3f} "
          f"(range {by_phylum['pearson_r_logsize'].min():.3f} to "
          f"{by_phylum['pearson_r_logsize'].max():.3f})")
    print(f"  better model: "
          f"{by_phylum['better_model'].value_counts().to_dict()}")
    print(f"  phylum-mean size vs phylum-mean KO count: "
          f"r = {conf.iloc[0]['pearson_r']:.3f} (p = {conf.iloc[0]['pearson_p']:.2e})")

    print("[4b/5] Linear vs log(size) model comparison")
    comp_long, comp_wide = correlation_comparison(ko_sum, meta, big_phyla)
    pooled = comp_wide[comp_wide["scope"] == "POOLED"].iloc[0]
    n_log = (comp_wide["better_model"] == "log_size").sum()
    print(f"  pooled: linear r={pooled['pearson_r_linear']:.3f} "
          f"(R2={pooled['ols_r2_linear']:.3f}) vs log(size) "
          f"r={pooled['pearson_r_log_size']:.3f} "
          f"(R2={pooled['ols_r2_log_size']:.3f}); "
          f"dAIC={pooled['dAIC_linear_minus_log']:.1f}")
    print(f"  better model is log_size in {n_log}/{len(comp_wide) - 1} phyla "
          f"(+ pooled)")

    print("[5/5] Writing outputs")
    ko_corr.round(5).to_csv(f"{WORK_DIR}/ko_size_correlation.tsv", sep="\t")
    summary.round(5).to_csv(f"{WORK_DIR}/sum_ko_size_correlation.tsv",
                            sep="\t", index=False)
    by_phylum.round(5).to_csv(f"{WORK_DIR}/ko_size_correlation_by_phylum.tsv",
                              sep="\t")
    conf.round(5).to_csv(f"{WORK_DIR}/phylum_mean_confounding.tsv",
                         sep="\t", index=False)
    comp_long.round(5).to_csv(
        f"{WORK_DIR}/correlation_linear_vs_log_long.tsv", sep="\t",
        index=False)
    comp_wide.round(5).to_csv(
        f"{WORK_DIR}/correlation_linear_vs_log_wide.tsv", sep="\t",
        index=False)

    for name, fig in [("fig1_sum_vs_size", fig1_regression(ko_sum, meta, summary)),
                      ("fig2_per_ko_correlations", fig2_per_ko(ko_corr)),
                      ("fig3_linear_vs_log_examples",
                       fig3_linear_vs_log(ko_sum, meta, list(big_phyla)))]:
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
