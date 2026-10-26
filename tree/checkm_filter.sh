#!/usr/bin/env bash
# Run CheckM (v1, lineage-specific markers) on a folder of genomes and keep those that meet:
#   completeness >= 90%, contamination <= 5%, genome size >= 0.1 Mb (100,000 bp)
#
# Usage:
#   ./checkm_filter.sh -g genome_dir [-o outdir] [-x fna] [-t threads] [-p pplacer_threads]
#       [-r] [-c 90] [-m 5] [-s 100000]
#
#   -g  folder of genome FASTA files
#   -o  output folder (default: checkm_out)
#   -x  file extension (default: fna)
#   -t  threads (default: all)
#   -p  pplacer threads (default: 4; each extra thread adds a lot of RAM)
#   -r  use --reduced_tree (~14 GB RAM instead of ~40 GB; slightly less precise lineage)
#   -c / -m / -s  thresholds: min completeness, max contamination, min size in bp
#
# Outputs (in outdir):
#   quality_report.tsv     CheckM extended QA table (incl. genome size)
#   quality_filtered.tsv   all genomes + PASS/FAIL + reason
#   pass_ids.txt           IDs of genomes that pass
#   pass_genomes/          hard links (or copies) of passing genomes, ready for GTDB-Tk
set -euo pipefail

GENOME_DIR=""
OUTDIR="checkm_out"
EXT="fna"
THREADS="$(nproc)"
PPLACER=4
REDUCED=()

MIN_COMP=90
MAX_CONT=5
MIN_SIZE=100000

usage() { sed -n '2,23p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
while getopts "g:o:x:t:p:rc:m:s:h" opt; do
  case "$opt" in
    g) GENOME_DIR="$OPTARG" ;; o) OUTDIR="$OPTARG" ;; x) EXT="${OPTARG#.}" ;;
    t) THREADS="$OPTARG" ;;    p) PPLACER="$OPTARG" ;; r) REDUCED=(--reduced_tree) ;;
    c) MIN_COMP="$OPTARG" ;;   m) MAX_CONT="$OPTARG" ;; s) MIN_SIZE="$OPTARG" ;; *) usage ;;
  esac
done
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: -g genome folder required" >&2; usage; }
command -v checkm >/dev/null || { echo "ERROR: checkm not in PATH (activate its conda env?)" >&2; exit 1; }
GENOME_DIR="$(cd "$GENOME_DIR" && pwd)"
mkdir -p "$OUTDIR"; OUTDIR="$(cd "$OUTDIR" && pwd)"

N=$(find -L "$GENOME_DIR" -maxdepth 1 -type f -name "*.$EXT" | wc -l)
echo ">> $N genomes (*.$EXT) in $GENOME_DIR"
[[ "$N" -gt 0 ]] || { echo "ERROR: no *.$EXT files found (broken symlinks? wrong -x?)" >&2; exit 1; }

# --- Run CheckM lineage workflow (skipped if it already finished) ---------------------
LW="$OUTDIR/lineage_wf"
if [[ ! -f "$OUTDIR/.lineage_wf.done" ]]; then
  rm -rf "$LW"   # CheckM cannot resume a partial run
  echo ">> Running checkm lineage_wf ($THREADS threads, $PPLACER pplacer threads) - may take hours"
  checkm lineage_wf -x "$EXT" -t "$THREADS" --pplacer_threads "$PPLACER" "${REDUCED[@]}" \
    --tab_table -f "$OUTDIR/lineage_wf_summary.tsv" "$GENOME_DIR" "$LW"
  touch "$OUTDIR/.lineage_wf.done"
fi

# --- Extended QA table (adds 'Genome size (bp)', N50, contig counts, GC, ...) ----------
REPORT="$OUTDIR/quality_report.tsv"
if [[ ! -s "$REPORT" ]]; then
  echo ">> Running checkm qa -o 2"
  checkm qa -o 2 --tab_table -t "$THREADS" -f "$REPORT" "$LW/lineage.ms" "$LW"
fi

# --- Filter (columns located by header name) ------------------------------------------
awk -F'\t' -v OFS='\t' -v c="$MIN_COMP" -v m="$MAX_CONT" -v s="$MIN_SIZE" \
    -v passf="$OUTDIR/pass_ids.txt" '
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i
          if (!("Bin Id" in h) || !("Completeness" in h) || !("Contamination" in h) || !("Genome size (bp)" in h)) {
            print "ERROR: unexpected CheckM header" > "/dev/stderr"; exit 1 }
          print "Name","Marker_lineage","Completeness","Contamination","Genome_Size_bp","Status","Reason"; next }
  { r=""
    if ($h["Completeness"]+0     <  c) r=r "completeness<" c ";"
    if ($h["Contamination"]+0    >  m) r=r "contamination>" m ";"
    if ($h["Genome size (bp)"]+0 <  s) r=r "size<" s "bp;"
    st = (r=="") ? "PASS" : "FAIL"
    if (st=="PASS") print $h["Bin Id"] > passf
    print $h["Bin Id"], $h["Marker lineage"], $h["Completeness"], $h["Contamination"], \
          $h["Genome size (bp)"], st, (r==""?"-":r) }
' "$REPORT" > "$OUTDIR/quality_filtered.tsv"
touch "$OUTDIR/pass_ids.txt"

# --- Collect passing genomes as real files (hard link, fall back to copy) -------------
PASS_DIR="$OUTDIR/pass_genomes"
rm -rf "$PASS_DIR"; mkdir -p "$PASS_DIR"
while read -r id; do
  real="$(readlink -f "$GENOME_DIR/$id.$EXT")"
  ln "$real" "$PASS_DIR/$id.$EXT" 2>/dev/null || cp "$real" "$PASS_DIR/$id.$EXT"
done < "$OUTDIR/pass_ids.txt"

NPASS=$(wc -l < "$OUTDIR/pass_ids.txt")
NTOT=$(( $(wc -l < "$REPORT") - 1 ))
echo ">> $NPASS / $NTOT genomes pass (completeness>=$MIN_COMP, contamination<=$MAX_CONT, size>=$MIN_SIZE bp)"
echo ">> Fail reasons:"
awk -F'\t' 'NR>1 && $6=="FAIL" { n=split($7,a,";"); for(i=1;i<n;i++) cnt[a[i]]++ }
            END { for (k in cnt) printf "     %-22s %d\n", k, cnt[k] }' "$OUTDIR/quality_filtered.tsv"
echo ">> Passing genomes: $PASS_DIR"
