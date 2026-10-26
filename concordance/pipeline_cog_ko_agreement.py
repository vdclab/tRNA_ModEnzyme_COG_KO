#!/usr/bin/env python3
"""
Pipeline: COG-KO agreement analysis for tRNA modification systems.

Inputs (tab-delimited, in INPUT_DIR):
  COG.txt  binary matrix, rows = 1322 genomes, cols = COG keys (1 present, 0 absent)
  KO.txt   binary matrix, rows = 1322 genomes, cols = KO keys
  key.txt  COG-KO pair mapping (header: COG, KO); some keys map to multiple partners

Outputs (in OUTPUT_DIR):
  agreement_matrix.tsv            genomes x unique COG-KO pairs (1 = agree, 0 = disagree)
  pair_agreement_stats.tsv        per-pair contingency counts and agreement %
  summary_by_mapping_category.tsv agreement by COG/KO multi-mapping category
  genome_phylum_map.tsv           genome -> NCBI taxid -> phylum (fetched from NCBI, cached)
  table_agreement_by_phylum.tsv   mean agreement per pair per phylum
  table_genome_outliers.tsv       bottom/top 2.5% genomes by mean agreement
  table_pca_top_loadings.tsv      pairs driving PC1/PC2
  table_phi_correlation.tsv       phi coefficient vs raw agreement per pair
  fig1_pair_agreement_bars.*      per-pair agreement, sorted
  fig2_disagreement_composition.* COG-only vs KO-only disagreement decomposition
  fig3_agreement_heatmap.*        genome x pair agreement heatmap
  fig4_agreement_by_phylum_heatmap.*  pair x phylum mean agreement
  fig5_genome_agreement_sorted.*  genomes sorted by agreement + phylum ribbon
  fig6_agreement_by_phylum_box.*  per-genome agreement by phylum
  fig7_pca_agreement.*            PCA of agreement matrix, colored by phylum
  fig8_phi_vs_agreement.*         phi correlation vs raw agreement

Usage: python3 pipeline_cog_ko_agreement.py
"""

import os
import re
import time
import json
import shutil
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.decomposition import PCA
from scipy.stats import spearmanr

# ---------------------------------------------------------------- configuration
INPUT_DIR = "/mnt/user-uploads"
OUTPUT_DIR = "/mnt/results/cog_ko_agreement"
WORK_DIR = "/workspace/cog_ko/pipeline_out"   # local scratch (S3 mount cannot overwrite)
MIN_PHYLUM_GENOMES = 10                        # phyla below this pooled as "Other (<10)"

plt.rcParams["font.family"] = ["Liberation Sans", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------- 1. load inputs
def load_inputs():
    cog = pd.read_csv(f"{INPUT_DIR}/COG.txt", sep="\t", index_col=0)
    ko = pd.read_csv(f"{INPUT_DIR}/KO.txt", sep="\t", index_col=0)
    key = pd.read_csv(f"{INPUT_DIR}/key.txt", sep="\t")

    common = cog.index.intersection(ko.index)
    cog, ko = cog.loc[common], ko.loc[common]
    log(f"COG matrix {cog.shape}, KO matrix {ko.shape}, {len(common)} shared genomes")

    # pairs without a KO cannot be compared
    no_ko = key[key["KO"].isna()]["COG"].tolist()
    key = key[key["KO"].notna()]
    # deduplicate COG-KO combinations (same pair listed for several reaction rows)
    n0 = len(key)
    key = key.drop_duplicates(subset=["COG", "KO"], keep="first").reset_index(drop=True)
    log(f"excluded {len(no_ko)} pairs without KO {no_ko}; "
        f"removed {n0 - len(key)} duplicate combinations; {len(key)} unique pairs")
    return cog, ko, key


# ---------------------------------------------------------------- 2. agreement matrix
def build_agreement(cog, ko, key):
    agree = pd.DataFrame(index=cog.index)
    records = []
    for _, row in key.iterrows():
        c, k = row["COG"], row["KO"]
        label = f"{c}_{k}"
        cv, kv = cog[c].astype(int), ko[k].astype(int)
        agree[label] = (cv == kv).astype(int)
        n11 = int(((cv == 1) & (kv == 1)).sum())
        n00 = int(((cv == 0) & (kv == 0)).sum())
        n10 = int(((cv == 1) & (kv == 0)).sum())   # COG present, KO absent
        n01 = int(((cv == 0) & (kv == 1)).sum())   # KO present, COG absent
        records.append(dict(pair=label, COG=c, KO=k, n_both_present=n11,
                            n_both_absent=n00, n_COG_only=n10, n_KO_only=n01,
                            n_agree=n11 + n00, n_disagree=n10 + n01,
                            agreement_pct=100 * (n11 + n00) / len(cv)))
    stats = pd.DataFrame(records)

    cog_counts = key["COG"].value_counts()
    ko_counts = key["KO"].value_counts()
    stats["COG_multi_mapped"] = stats["COG"].map(cog_counts) > 1
    stats["KO_multi_mapped"] = stats["KO"].map(ko_counts) > 1
    return agree, stats


# ---------------------------------------------------------------- 3. phylum mapping
def phylum_from_taxon_xml(xml_text):
    root = ET.fromstring(xml_text)
    out = {}
    for tr in root.findall(".//Taxon"):
        tid = tr.findtext("TaxId")
        p = "unresolved"
        for node in tr.findall(".//LineageEx/Taxon"):
            if node.findtext("Rank") in ("phylum", "superphylum"):
                p = node.findtext("ScientificName")
                break
        out[tid] = p
    return out


def eutils(url, params=None, data=None, retries=5):
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, data=data, timeout=180) as r:
                return r.read().decode()
        except Exception as e:
            log(f"  retry {attempt}: {str(e)[:80]}")
            time.sleep(15)
    raise RuntimeError(f"eutils failed: {url}")


def fetch_phyla(genome_ids, cache_path):
    # cached result from a previous run avoids re-hitting the NCBI API
    if os.path.exists(cache_path):
        ph = pd.read_csv(cache_path, sep="\t", index_col=0)["phylum"]
        if set(ph.index) >= set(genome_ids):
            log(f"phylum map loaded from cache ({len(ph)} genomes)")
            return ph.loc[genome_ids]

    accs = list(genome_ids)
    acc2uid = {}
    for i in range(0, len(accs), 100):
        b = accs[i:i + 100]
        term = " OR ".join(f"{a}[Assembly Accession]" for a in b)
        res = json.loads(eutils("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                                data=urllib.parse.urlencode(
                                    {"db": "assembly", "term": term,
                                     "retmax": 200, "retmode": "json"}).encode()))
        uids = res["esearchresult"]["idlist"]
        if uids:
            summ = json.loads(eutils(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                {"db": "assembly", "id": ",".join(uids), "retmode": "json"}))["result"]
            for uid in summ.get("uids", []):
                m = re.match(r"GC[AF]_(\d+)", summ[uid].get("assemblyaccession") or "")
                if m:
                    acc2uid[m.group(1)] = uid
        time.sleep(1)
    log(f"assemblies resolved: {len(acc2uid)}")

    taxids = {}
    uid_list = [acc2uid[a.split("_")[1].split(".")[0]] for a in accs
                if a.split("_")[1].split(".")[0] in acc2uid]
    for i in range(0, len(uid_list), 200):
        b = uid_list[i:i + 200]
        summ = json.loads(eutils(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            {"db": "assembly", "id": ",".join(b), "retmode": "json"}))["result"]
        for uid in summ.get("uids", []):
            taxids[uid] = str(summ[uid].get("taxid"))
        time.sleep(1)

    phylum = {}
    tids = list(set(taxids.values()))
    for i in range(0, len(tids), 150):
        b = tids[i:i + 150]
        xml = eutils("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                     {"db": "taxonomy", "id": ",".join(b), "retmode": "xml"})
        phylum.update(phylum_from_taxon_xml(xml))
        time.sleep(1)

    rows = [(a, taxids.get(acc2uid.get(a.split("_")[1].split(".")[0]), ""),
             phylum.get(taxids.get(acc2uid.get(a.split("_")[1].split(".")[0]), ""), "unresolved"))
            for a in accs]
    ph = pd.DataFrame(rows, columns=["genome", "taxid", "phylum"]).set_index("genome")["phylum"]
    ph.to_csv(cache_path, sep="\t")
    log("phylum distribution:\n" + ph.value_counts().head(10).to_string())
    return ph


# ---------------------------------------------------------------- 4. figures
def fig1_bars(stats, out):
    s = stats.sort_values("agreement_pct", ascending=False).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(s))
    ax.bar(x, s["agreement_pct"], color="#0279EE", width=0.8)
    ax.axhline(90, color="#FF9400", ls="--", lw=1, label="90%")
    ax.axhline(stats["agreement_pct"].mean(), color="black", ls=":", lw=1,
               label=f"mean = {stats['agreement_pct'].mean():.1f}%")
    ax.set_xticks(x[::5])
    ax.set_xticklabels(s["pair"][::5], rotation=90, fontsize=7)
    ax.set_ylabel("Genomes with COG–KO agreement (%)")
    ax.set_ylim(0, 102)
    ax.set_xlabel("COG–KO pair (sorted by agreement)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{out}/fig1_pair_agreement_bars.{ext}",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)


def fig2_disagreement(stats, out):
    s = stats.sort_values("n_disagree", ascending=False).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(s))
    ax.bar(x, s["n_COG_only"], color="#FF9400", label="COG present, KO absent")
    ax.bar(x, s["n_KO_only"], bottom=s["n_COG_only"], color="#75A025",
           label="KO present, COG absent")
    ax.set_xticks(x[::5])
    ax.set_xticklabels(s["pair"][::5], rotation=90, fontsize=7)
    ax.set_ylabel("Genomes with disagreement (n)")
    ax.set_xlabel("COG–KO pair (sorted by disagreement count)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{out}/fig2_disagreement_composition.{ext}",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)


def fig3_heatmap(agree, out):
    fig, ax = plt.subplots(figsize=(16, 10))
    im = ax.imshow(agree.values, aspect="auto", interpolation="nearest",
                   cmap=ListedColormap(["#ECE9E2", "#0279EE"]))
    ax.set_xticks(range(agree.shape[1]))
    ax.set_xticklabels(agree.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(0, agree.shape[0], 100))
    ax.set_yticklabels(agree.index[::100], fontsize=6)
    ax.set_xlabel("COG–KO pair")
    ax.set_ylabel("Genome")
    ax.set_title("COG–KO agreement per genome (blue = agree)")
    fig.colorbar(im, ax=ax, shrink=0.5, ticks=[0.25, 0.75]).set_ticklabels(
        ["disagree", "agree"])
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{out}/fig3_agreement_heatmap.{ext}",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)


def fig4_phylum_heatmap(agree, ph_maj, out):
    rows = []
    for pair in agree.columns:
        row = {"pair": pair}
        for p, v in agree[pair].groupby(ph_maj):
            row[p] = 100 * v.mean()
        rows.append(row)
    by_ph = pd.DataFrame(rows).set_index("pair")
    order = [p for p in by_ph.columns if p != "Other (<10)"] + ["Other (<10)"]
    by_ph = by_ph[order].sort_values("Pseudomonadota", ascending=False)

    fig, ax = plt.subplots(figsize=(9, 14))
    im = ax.imshow(by_ph.values, aspect="auto", cmap="RdYlGn", vmin=0, vmax=100)
    ax.set_xticks(range(by_ph.shape[1]))
    ax.set_xticklabels(by_ph.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(by_ph.shape[0]))
    ax.set_yticklabels(by_ph.index, fontsize=5)
    ax.set_title("Mean COG–KO agreement (%) by phylum")
    fig.colorbar(im, ax=ax, shrink=0.4, label="% genomes agreeing")
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{out}/fig4_agreement_by_phylum_heatmap.{ext}",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)
    by_ph.round(1).to_csv(f"{out}/table_agreement_by_phylum.tsv", sep="\t")


def fig5_6_genome_level(agree, ph_maj, out):
    genome_agree = agree.mean(axis=1) * 100
    gs = pd.DataFrame({"mean_agreement": genome_agree, "phylum": ph_maj})
    gs_sorted = gs.sort_values("mean_agreement")

    lo, hi = gs_sorted["mean_agreement"].quantile([0.025, 0.975])
    gs_sorted[(gs_sorted["mean_agreement"] <= lo) |
              (gs_sorted["mean_agreement"] >= hi)].to_csv(
        f"{out}/table_genome_outliers.tsv", sep="\t")

    cmap = plt.get_cmap("tab20")
    phyla = sorted(gs["phylum"].unique())
    colors = {p: cmap(i / max(1, len(phyla) - 1)) for i, p in enumerate(phyla)}

    # fig5: sorted bars + phylum ribbon
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(16, 9), height_ratios=[4, 0.35],
                                  sharex=True, constrained_layout=True)
    x = np.arange(len(gs_sorted))
    ax.bar(x, gs_sorted["mean_agreement"],
           color=[colors[p] for p in gs_sorted["phylum"]], width=1.0)
    ax.axhline(gs["mean_agreement"].mean(), color="black", ls=":", lw=1,
               label=f"mean = {gs['mean_agreement'].mean():.1f}%")
    ax.set_ylabel("Mean COG–KO agreement per genome (%)")
    ax.set_xlabel("Genomes sorted by mean agreement")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="lower right")
    for i, p in enumerate(gs_sorted["phylum"]):
        ax2.bar(i, 1, color=colors[p], width=1.0)
    ax2.set_ylabel("Phylum")
    ax2.set_yticks([])
    for ext in ("svg", "png"):
        fig.savefig(f"{out}/fig5_genome_agreement_sorted.{ext}",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)

    # fig6: boxplot by phylum
    fig, ax = plt.subplots(figsize=(10, 5))
    data = [gs.loc[gs["phylum"] == p, "mean_agreement"].values for p in phyla]
    bp = ax.boxplot(data, tick_labels=phyla, patch_artist=True, showfliers=False)
    for patch, p in zip(bp["boxes"], phyla):
        patch.set_facecolor(colors[p])
    rng = np.random.default_rng(1)
    for i, p in enumerate(phyla):
        y = gs.loc[gs["phylum"] == p, "mean_agreement"]
        ax.scatter(rng.normal(i + 1, 0.06, len(y)), y, s=4, color="black",
                   alpha=0.25, zorder=3)
    ax.set_ylabel("Mean COG–KO agreement per genome (%)")
    ax.tick_params(axis="x", rotation=45)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{out}/fig6_agreement_by_phylum_box.{ext}",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)
    return gs


def fig7_pca(agree, gs, out):
    pca = PCA(n_components=3)
    emb = pca.fit_transform(agree.values)
    evr = pca.explained_variance_ratio_ * 100

    cmap = plt.get_cmap("tab20")
    phyla = sorted(gs["phylum"].unique())
    colors = {p: cmap(i / max(1, len(phyla) - 1)) for i, p in enumerate(phyla)}

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax, (i, j) in zip(axes, [(0, 1), (0, 2)]):
        for p in phyla:
            m = (gs["phylum"] == p).values
            ax.scatter(emb[m, i], emb[m, j], s=8, alpha=0.5, color=colors[p], label=p)
        ax.set_xlabel(f"PC{i+1} ({evr[i]:.1f}%)")
        ax.set_ylabel(f"PC{j+1} ({evr[j]:.1f}%)")
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=10,
               loc="center left", bbox_to_anchor=(1.01, 0.5),
               title="Phylum", title_fontsize=11)
    fig.suptitle("PCA of COG–KO agreement profiles")
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{out}/fig7_pca_agreement.{ext}", bbox_inches="tight",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)

    load = pd.DataFrame({"pair": agree.columns,
                         "PC1": pca.components_[0], "PC2": pca.components_[1]})
    load["absPC1"] = load["PC1"].abs()
    load.sort_values("absPC1", ascending=False).head(15).to_csv(
        f"{out}/table_pca_top_loadings.tsv", sep="\t", index=False)


def fig8_phi(cog, ko, key, stats, out):
    def phi(a, b):
        a, b = np.asarray(a, float), np.asarray(b, float)
        if a.std() == 0 or b.std() == 0:
            return np.nan
        return np.corrcoef(a, b)[0, 1]

    phi_rows = [dict(pair=f"{r['COG']}_{r['KO']}", phi=phi(cog[r["COG"]], ko[r["KO"]]))
                for _, r in key.iterrows()]
    phi_df = pd.DataFrame(phi_rows).set_index("pair")
    phi_df = phi_df[~phi_df.index.duplicated(keep="first")]
    stats_u = stats.drop_duplicates(subset="pair", keep="first").set_index("pair")
    phi_df["agreement_pct"] = stats_u["agreement_pct"]
    phi_df.sort_values("phi").to_csv(f"{out}/table_phi_correlation.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(phi_df["phi"], phi_df["agreement_pct"], s=12, color="#0279EE", alpha=0.6)
    for p in ["COG0037_K21947", "COG0564_K06175", "COG0042_K05541", "COG0566_K00556",
              "COG0301_K03151", "COG0336_K00554"]:
        if p in phi_df.index:
            ax.annotate(p, (phi_df.loc[p, "phi"], phi_df.loc[p, "agreement_pct"]),
                        fontsize=6, alpha=0.8)
    ax.set_xlabel("Phi correlation (COG vs KO presence)")
    ax.set_ylabel("Agreement (%)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{out}/fig8_phi_vs_agreement.{ext}",
                    dpi=200 if ext == "png" else None)
    plt.close(fig)


# ---------------------------------------------------------------- 5. main
def main():
    os.makedirs(WORK_DIR, exist_ok=True)
    cog, ko, key = load_inputs()
    agree, stats = build_agreement(cog, ko, key)

    agree.to_csv(f"{WORK_DIR}/agreement_matrix.tsv", sep="\t")
    stats.to_csv(f"{WORK_DIR}/pair_agreement_stats.tsv", sep="\t", index=False)
    summary = stats.groupby(["COG_multi_mapped", "KO_multi_mapped"]).agg(
        n_pairs=("pair", "count"),
        mean_agreement_pct=("agreement_pct", "mean"),
        median_agreement_pct=("agreement_pct", "median")).round(1)
    summary.to_csv(f"{WORK_DIR}/summary_by_mapping_category.tsv", sep="\t")
    log(f"agreement matrix {agree.shape}; mean agreement "
        f"{stats['agreement_pct'].mean():.1f}%")

    ph = fetch_phyla(list(agree.index), f"{WORK_DIR}/genome_phylum_map.tsv")
    counts = ph.value_counts()
    major = counts[counts >= MIN_PHYLUM_GENOMES].index.tolist()
    ph_maj = ph.where(ph.isin(major), "Other (<10)")

    fig1_bars(stats, WORK_DIR)
    fig2_disagreement(stats, WORK_DIR)
    fig3_heatmap(agree, WORK_DIR)
    fig4_phylum_heatmap(agree, ph_maj, WORK_DIR)
    gs = fig5_6_genome_level(agree, ph_maj, WORK_DIR)
    fig7_pca(agree, gs, WORK_DIR)
    fig8_phi(cog, ko, key, stats, WORK_DIR)

    # sync to the results mount (remove-then-copy: the mount forbids overwrite)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for f in os.listdir(WORK_DIR):
        dst = os.path.join(OUTPUT_DIR, f)
        if os.path.exists(dst):
            os.remove(dst)
        shutil.copy(os.path.join(WORK_DIR, f), dst)
    log(f"done — all outputs in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
