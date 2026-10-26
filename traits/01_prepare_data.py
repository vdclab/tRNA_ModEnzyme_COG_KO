#!/usr/bin/env python3
"""Step 1 (v2 robustness rerun): data preparation with updated inputs.

Changes vs v1:
  - select_trait.txt / pathway_groups.txt / genome_traits.txt from updated
    uploads (Phylum and Habitats dropped from trait list; Oxygen_Req. now
    Categorical; pathway_groups has NO within_pathway_independent column)
  - tree pre-rooted at Aquificota (bac120_aquificota_rooted.tree)
  - outputs to /mnt/results/ko_trait_analysis_v2/00_processed/
  - Phylum kept as metadata column for the per-phylum sensitivity analysis
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

UP = Path("/mnt/user-uploads")
IN = UP / "input"
OUT = Path("/mnt/results/ko_trait_analysis_v2/00_processed")
OUT.mkdir(parents=True, exist_ok=True)

ko = pd.read_csv(IN / "KO_AP.txt", sep="\t", index_col=0)
traits = pd.read_csv(UP / "genome_traits.txt", sep="\t", index_col=0, dtype=str)
tree_newick = Path("/workspace/ko_trait_analysis_v2/bac120_aquificota_rooted.tree"
                   ).read_text().strip()

sel = pd.read_csv(UP / "select_trait.txt", sep="\t", header=None,
                  names=["trait", "type"])
pg = pd.read_csv(UP / "pathway_groups.txt", sep="\t")
pg.columns = [c.strip() for c in pg.columns]

checks = []
checks.append(f"KO matrix: {ko.shape[0]} genomes x {ko.shape[1]} KOs; "
              f"values = {sorted(pd.unique(ko.values.ravel()))}")
checks.append(f"Traits table: {traits.shape[0]} genomes x {traits.shape[1]} columns")

tips = set(re.findall(r"[(,]([A-Za-z0-9_.]+):", tree_newick))
assert tips == set(ko.index) == set(traits.index), "genome ID mismatch across files"
checks.append("ID consistency: rooted tree tips == KO matrix rows == traits rows "
              f"({len(tips)})")

missing_kos = set(ko.columns) - set(pg["KO"])
extra_kos = set(pg["KO"]) - set(ko.columns)
checks.append(f"pathway_groups KOs not in matrix: {sorted(missing_kos)}")
checks.append(f"matrix KOs not in pathway_groups: {sorted(extra_kos)}")
checks.append(f"within_pathway_independent column present: "
              f"{'within_pathway_independent' in pg.columns}")

pg = pg[pg["KO"].isin(ko.columns)].copy()
pathway_members = pg.groupby("pathway")["KO"].apply(list)
completeness = pd.DataFrame(index=ko.index)
for pathway, kos in pathway_members.items():
    completeness[pathway] = ko[kos].mean(axis=1)

zero_var = [p for p in completeness.columns if completeness[p].nunique() == 1]
checks.append(f"Pathways: {completeness.shape[1]}; zero-variance pathways: {zero_var}")

# trait table: selected traits + Phylum as metadata (for per-phylum analysis)
trait_tbl = pd.DataFrame(index=traits.index)
trait_tbl["Phylum"] = traits["Phylum"].replace("", np.nan)
for _, r in sel.iterrows():
    t = r["trait"]
    if t not in traits.columns:
        checks.append(f"WARNING: trait '{t}' not found in genome_traits.txt")
        continue
    trait_tbl[t] = traits[t].replace("", np.nan)

binary_maps = {
    "Gram_Stain": {"positive": 1, "negative": 0},
    "Motility": {"motile": 1, "nonmotile": 0},
}
for t, m in binary_maps.items():
    if t in trait_tbl.columns:
        unmapped = set(trait_tbl[t].dropna()) - set(m)
        if unmapped:
            checks.append(f"WARNING: unmapped levels for {t}: {unmapped}")
        trait_tbl[t + "_coded"] = trait_tbl[t].map(m)

for t in sel["trait"]:
    if t in trait_tbl.columns:
        n = trait_tbl[t].notna().sum()
        checks.append(f"Trait {t}: type={sel.loc[sel.trait == t, 'type'].iloc[0]}, "
                      f"complete cases={n}, levels={trait_tbl[t].nunique()}")

ko.to_csv(OUT / "ko_binary.csv")
completeness.to_csv(OUT / "pathway_completeness.csv")
trait_tbl.to_csv(OUT / "traits_selected.csv")
pathway_members.rename("KOs").to_csv(OUT / "pathway_members.csv")
(OUT / "summary_checks.txt").write_text("\n".join(checks) + "\n")

print("\n".join(checks))
print(f"\nCompleteness matrix: {completeness.shape}")
