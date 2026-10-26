#!/usr/bin/env bash
# Workflow: copy-number distribution of COGs and KOs, overall and by phylum.
#
# Usage:
#   bash run_count_distribution.sh [COG_count.txt] [KO_count.txt] [genome_phylum.txt] [outdir] [min_genomes]
#
# Defaults: COG_count.txt KO_count.txt genome_phylum.txt results 10
set -euo pipefail

COG=${1:-COG_count.txt}
KO=${2:-KO_count.txt}
PHY=${3:-genome_phylum.txt}
OUT=${4:-results}
MIN=${5:-10}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SCRIPT_DIR}/count_distribution.py"

# ---- checks -----------------------------------------------------------------
for f in "$COG" "$KO" "$PHY" "$PY"; do
    [[ -s "$f" ]] || { echo "Error: missing or empty file: $f" >&2; exit 1; }
done
python3 -c "import pandas, numpy, matplotlib" 2>/dev/null || {
    echo "Error: python3 with pandas, numpy and matplotlib is required" >&2; exit 1; }

mkdir -p "$OUT/clean"

# ---- 0) normalise inputs: remove Windows line endings, ensure final newline --
clean() { sed 's/\r$//' "$1" | sed -e '$a\' > "$2"; }
clean "$COG" "$OUT/clean/COG_count.txt"
clean "$KO"  "$OUT/clean/KO_count.txt"
clean "$PHY" "$OUT/clean/genome_phylum.txt"

echo "Genomes: COG=$(($(wc -l < "$OUT/clean/COG_count.txt") - 1))," \
     "KO=$(($(wc -l < "$OUT/clean/KO_count.txt") - 1))," \
     "phylum file=$(($(wc -l < "$OUT/clean/genome_phylum.txt") - 1))"

# ---- 1-3) distributions, grouping by phylum, heatmaps -----------------------
python3 "$PY" "$OUT/clean/COG_count.txt" "$OUT/clean/genome_phylum.txt" COG "$OUT" "$MIN"
python3 "$PY" "$OUT/clean/KO_count.txt"  "$OUT/clean/genome_phylum.txt" KO  "$OUT" "$MIN"

rm -r "$OUT/clean"
echo "Done. Outputs in $OUT/:"
ls -1 "$OUT"
