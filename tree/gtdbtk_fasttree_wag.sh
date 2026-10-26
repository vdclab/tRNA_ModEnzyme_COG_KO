#!/bin/bash

# Build a phylogenomic tree of your genomes with GTDB-Tk:
#   identify (find 120 bacterial / 53 archaeal marker genes)
#   -> align (concatenated, trimmed protein MSA of YOUR genomes only)
#   -> infer (FastTree; optional IQ-TREE for a slower, more accurate ML tree with bootstraps)
#   -> optional classify_wf (GTDB taxonomy + iTOL label file for the tree)
# Each step is skipped if already completed, so you can rerun after an interruption.
#
# Usage:
#   ./gtdbtk_tree.sh -g genomes_dir [-o outdir] [-c cpus] [-x ext] [-m WAG|LG|JTT] [-C] [-I]
#
#   -g  folder containing genome FASTA files (one genome per file; .gz is fine)
#   -o  output directory (default: gtdbtk_out)
#   -c  CPUs (default: all)
#   -x  file extension, e.g. fna, fna.gz, fa, fasta (default: auto-detect)
#   -m  FastTree protein model for 'gtdbtk infer' (default: WAG)
#   -C  also run classify_wf to get GTDB taxonomy + iTOL labels (slow, ~65+ GB RAM)
#   -I  also build an IQ-TREE ML tree with 1000 ultrafast bootstraps (needs iqtree2)
set -euo pipefail

GENOME_DIR=$1
OUTDIR="gtdbtk_out"
CPUS=12
EXT=fna
MODEL="WAG"; CLASSIFY=0; IQTREE=0

usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
while getopts "g:o:c:x:m:CIh" opt; do
  case "$opt" in
    g) GENOME_DIR="$OPTARG" ;; o) OUTDIR="$OPTARG" ;; c) CPUS="$OPTARG" ;;
    x) EXT="${OPTARG#.}" ;;    m) MODEL="$OPTARG" ;;  C) CLASSIFY=1 ;; I) IQTREE=1 ;;
    *) usage ;;
  esac
done
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: -g genome directory required" >&2; usage; }
log() { echo "[$(date '+%F %T')] $*"; }

# --- Sanity checks ---------------------------------------------------------------------
source ~/data/bin/conda/etc/profile.d/conda.sh

conda activate gtdbtk-2.7.2

command -v gtdbtk >/dev/null || { echo "ERROR: gtdbtk not in PATH (activate its conda env?)" >&2; exit 1; }

[[ -n "${GTDBTK_DATA_PATH:-}" && -d "$GTDBTK_DATA_PATH" ]] || {
  echo "ERROR: GTDBTK_DATA_PATH is not set or not a directory (point it at the GTDB-Tk reference data)" >&2; exit 1; }
log "$(gtdbtk --version 2>&1 | grep -i 'gtdb-tk' | head -1 || true)"

if [[ -z "$EXT" ]]; then
  for e in fna.gz fa.gz fasta.gz fna fa fasta fas; do
    if compgen -G "$GENOME_DIR/*.$e" >/dev/null; then EXT="$e"; break; fi
  done
  [[ -n "$EXT" ]] || { echo "ERROR: no FASTA files found in $GENOME_DIR (use -x)" >&2; exit 1; }
fi
N=$(find -L "$GENOME_DIR" -maxdepth 1 -type f -name "*.$EXT" | wc -l)
log "Found $N genomes with extension .$EXT in $GENOME_DIR (using $CPUS CPUs)"
if find "$GENOME_DIR" -maxdepth 1 -name "*.$EXT" | grep -q '[[:space:]]'; then
  echo "ERROR: file names contain spaces; GTDB-Tk uses file names as genome IDs - rename them" >&2; exit 1
fi

mkdir -p "$OUTDIR/trees" "$OUTDIR/tmp"
OUTDIR="$(cd "$OUTDIR" && pwd)"
TMP="$OUTDIR/tmp"

# --- 1. identify -----------------------------------------------------------------------
if [[ ! -f "$OUTDIR/.identify.done" ]]; then
  log "Step 1/3: gtdbtk identify"
  gtdbtk identify --genome_dir "$GENOME_DIR" --extension "$EXT" \
    --out_dir "$OUTDIR/identify" --cpus "$CPUS" --tmpdir "$TMP"
  touch "$OUTDIR/.identify.done"
else log "Step 1/3: identify already done, skipping"; fi

# --- 2. align --------------------------------------------------------------------------
if [[ ! -f "$OUTDIR/.align.done" ]]; then
  log "Step 2/3: gtdbtk align"
  gtdbtk align --identify_dir "$OUTDIR/identify" --out_dir "$OUTDIR/align" \
    --cpus "$CPUS" --tmpdir "$TMP"
  touch "$OUTDIR/.align.done"
else log "Step 2/3: align already done, skipping"; fi

FILTERED=$(find "$OUTDIR/align" -name 'gtdbtk.*.filtered.tsv' -exec cat {} + 2>/dev/null | wc -l)
[[ "$FILTERED" -gt 0 ]] && log "NOTE: $FILTERED genomes dropped for too few marker AAs (see align/*filtered.tsv)"

# --- 3. infer (one tree per domain present: bac120 and/or ar53) ------------------------
BUILT=0
for dom in bac120; do
  msa_gz=$(find "$OUTDIR/align" -name "gtdbtk.${dom}.user_msa.fasta*" | head -1)
  [[ -z "$msa_gz" ]] && continue
  msa="$OUTDIR/trees/${dom}.user_msa.faa"
  if [[ "$msa_gz" == *.gz ]]; then gunzip -c "$msa_gz" > "$msa"; else cp "$msa_gz" "$msa"; fi
  ntaxa=$(grep -c '^>' "$msa")
  if [[ "$ntaxa" -lt 3 ]]; then log "Skipping $dom: only $ntaxa genomes"; continue; fi

  if [[ ! -f "$OUTDIR/trees/${dom}.fasttree.unrooted.tree" ]]; then
    log "Step 3/3: gtdbtk infer ($dom, $ntaxa genomes, FastTree $MODEL+GAMMA)"
    gtdbtk infer --msa_file "$msa" --out_dir "$OUTDIR/infer_$dom" \
      --prot_model "$MODEL" --gamma --cpus "$CPUS" --tmpdir "$TMP"
    tree=$(find "$OUTDIR/infer_$dom" -name '*unrooted.tree' | head -1)
    cp "$tree" "$OUTDIR/trees/${dom}.fasttree.unrooted.tree"
  fi
  BUILT=1

  if [[ "$IQTREE" -eq 1 ]]; then
    IQ=$(command -v iqtree2 || command -v iqtree || true)
    if [[ -z "$IQ" ]]; then log "WARNING: iqtree2 not found, skipping IQ-TREE"
    elif [[ ! -f "$OUTDIR/iqtree_$dom/$dom.treefile" ]]; then
      log "IQ-TREE ($dom): LG+F+G4, 1000 UFBoot (this can take many hours for ~1k genomes)"
      mkdir -p "$OUTDIR/iqtree_$dom"
      "$IQ" -s "$msa" -m LG+F+G4 -B 1000 -T "$CPUS" -pre "$OUTDIR/iqtree_$dom/$dom" -redo
      cp "$OUTDIR/iqtree_$dom/$dom.treefile" "$OUTDIR/trees/${dom}.iqtree.tree"
    fi
  fi
done
[[ "$BUILT" -eq 1 ]] || { echo "ERROR: no MSA produced - check $OUTDIR/identify logs" >&2; exit 1; }

<<'####'
# --- Optional: GTDB classification + iTOL label file ----------------------------------
if [[ "$CLASSIFY" -eq 1 ]]; then
  if [[ ! -f "$OUTDIR/.classify.done" ]]; then
    log "classify_wf (needs lots of RAM; pplacer step is the bottleneck)"
    EXTRA=()
    gtdbtk classify_wf -h 2>&1 | grep -q -- '--skip_ani_screen' && EXTRA+=(--skip_ani_screen)
    gtdbtk classify_wf --genome_dir "$GENOME_DIR" --extension "$EXT" \
      --out_dir "$OUTDIR/classify" --cpus "$CPUS" --pplacer_cpus "$(( CPUS < 8 ? CPUS : 8 ))" \
      --tmpdir "$TMP" "${EXTRA[@]}"
    touch "$OUTDIR/.classify.done"
  fi
  # iTOL label file: tip label = genome ID + lowest assigned GTDB rank
  LABELS="$OUTDIR/trees/itol_labels_gtdb.txt"
  { printf 'LABELS\nSEPARATOR TAB\nDATA\n'
    for f in "$OUTDIR"/classify/*summary.tsv "$OUTDIR"/classify/classify/*summary.tsv; do
      [[ -f "$f" ]] || continue
      awk -F'\t' 'NR>1 && $1!="user_genome" {
        n=split($2,r,";"); lab=$1
        for(i=n;i>=1;i--){ if(length(r[i])>3){ lab=$1" | "r[i]; break } }
        print $1"\t"lab }' "$f"
    done; } > "$LABELS"
  cp "$OUTDIR"/classify/*summary.tsv "$OUTDIR/trees/" 2>/dev/null || true
  log "GTDB taxonomy + iTOL labels written to $OUTDIR/trees/"
fi
####

rm -rf "$TMP"
log "Done. Trees are in $OUTDIR/trees/:"
ls -1 "$OUTDIR/trees/"


