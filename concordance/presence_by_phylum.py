#!/usr/bin/env python3
"""
Summarize presence of each COG and each KO by phylum.

Cell value (%) = (number of genomes in that phylum with the key present)
                 / (total number of genomes in the matrix) * 100

Usage:
    python presence_by_phylum.py COG.txt KO.txt genome_phylum.txt [outdir]

Outputs:
    COG_presence_summary.tsv   rows = COG keys, columns = phyla
    KO_presence_summary.tsv    rows = KO keys,  columns = phyla
"""
import os
import sys
import pandas as pd


def read_matrix(path):
    m = pd.read_csv(path, sep="\t", index_col=0, dtype=str)
    m.columns = m.columns.str.strip()
    m.index = m.index.astype(str).str.strip()
    return m.apply(pd.to_numeric, errors="raise").astype(int)


def read_phylum(path):
    p = pd.read_csv(path, sep="\t", dtype=str)
    p = p.apply(lambda s: s.str.strip())
    return p.set_index(p.columns[0])[p.columns[1]]


def summarize(matrix, phylum, key_name, out_path):
    ph = phylum.reindex(matrix.index).fillna("Unknown")
    n_missing = (ph == "Unknown").sum()
    if n_missing:
        print(f"Warning ({key_name}): {n_missing} genome(s) without phylum -> 'Unknown'",
              file=sys.stderr)

    total = len(matrix)                                   # e.g. 1322
    counts = matrix.groupby(ph).sum().T                   # keys x phyla
    order = ph.value_counts().index                       # phyla, largest first
    summary = counts[order] / total * 100                 # percentage
    summary.index.name = key_name
    summary.columns.name = None
    summary.to_csv(out_path, sep="\t")
    print(f"{key_name}: {summary.shape[0]} keys x {summary.shape[1]} phyla, "
          f"denominator = {total} genomes -> {out_path}", file=sys.stderr)


def main(cog_file, ko_file, phylum_file, outdir="."):
    os.makedirs(outdir, exist_ok=True)
    phylum = read_phylum(phylum_file)
    summarize(read_matrix(cog_file), phylum, "COG",
              os.path.join(outdir, "COG_presence_summary.tsv"))
    summarize(read_matrix(ko_file), phylum, "KO",
              os.path.join(outdir, "KO_presence_summary.tsv"))


if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    main(*sys.argv[1:])
