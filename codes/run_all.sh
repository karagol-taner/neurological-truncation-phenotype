#!/bin/bash
set -euo pipefail
C=$(cd "$(dirname "$0")" && pwd)
ROOT=$(dirname "$C")
[ $# -eq 2 ] || { echo "usage: bash codes/run_all.sh WORK_DIR LCR_BED_GZ"; exit 2; }
mkdir -p "$1"
W=$(cd "$1" && pwd)
BED=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
for f in variant_summary.txt.gz submission_summary.txt.gz gene_condition_source_id; do
  [ -s "$ROOT/$f" ] || { echo "missing $ROOT/$f"; exit 2; }
done
[ -s "$BED" ] || { echo "missing $BED"; exit 2; }
mkdir -p "$W"/{gate,pull,clinvar,run,out,gw/codes}
stamp() { echo "$(date -u +%H:%M:%S) $*"; }

python3 "$C/gate.py" "$W/gate" --jobs 2
mv "$C/gate_brain.json" "$W/out/gate_brain.json"
python3 "$C/grouping_report.py" "$W/gate" "$W/out/grouping_names.json"
stamp "gate done"

python3 "$C/panel_symbols.py" "$W/gate" "$W/clinvar/panel_symbols.json"
python3 "$C/vs_filter.py" "$W/clinvar/panel_symbols.json" "$ROOT/variant_summary.txt.gz" "$W/clinvar/variant_summary.txt.gz"
python3 "$C/ss_filter.py" "$W/clinvar/variant_summary.txt.gz" "$ROOT/submission_summary.txt.gz" "$W/clinvar/submission_summary.txt.gz"
stamp "ClinVar extract done"

python3 "$C/gnomad_union.py" "$W/gate" "$W/pull"
cp "$W/gate/strict/brain_whole/codes/derived_v2.json" "$W/pull/derived_v2.json"
cp "$C/gnomad_pull.py" "$W/pull/gnomad_pull.py"
python3 "$C/gnomad_pull_durable.py" "$W/pull"
stamp "gnomAD pull done"

SYS="$C/panel_v2.py"
setup_run() {
  local d="$W/run/$1"
  rm -rf "$d"; mkdir -p "$d/codes"
  cp "$C"/*.py "$C"/*.sh "$C/uniprot_rows.json" "$d/codes/"
  cp "$2" "$d/codes/panel_v2.py"
  cp "$W/pull/gnomad_plof.json" "$d/codes/"
  if [ -s "$W/tx_cache.json" ]; then cp "$W/tx_cache.json" "$d/codes/"; fi
  ln -s "$W/clinvar/variant_summary.txt.gz" "$d/variant_summary.txt.gz"
  ln -s "$W/clinvar/submission_summary.txt.gz" "$d/submission_summary.txt.gz"
  ln -s "$ROOT/gene_condition_source_id" "$d/gene_condition_source_id"
}
run() { bash "$W/run/$1/codes/run_unit.sh" "$W/run/$1" "$2" "$3" "$SYS"; stamp "$1 done"; }

setup_run brain_whole "$C/panel_v2.py"
run brain_whole 10 routes
[ -s "$W/tx_cache.json" ] || cp "$W/run/brain_whole/codes/tx_cache.json" "$W/tx_cache.json"
for u in epilepsy neurodevelopmental neurodegeneration malformation; do
  setup_run "brain_$u" "$ROOT/units/brain_$u/codes/panel_v2.py"
  run "brain_$u" 10 routes
done
setup_run sens_whole_kept "$C/panel_v2.py"
run sens_whole_kept 0 names
setup_run sens_whole_names10 "$C/panel_v2.py"
run sens_whole_names10 10 names

bash "$C/run_restrictions.sh" "$W/run" "$BED" "$W/clinvar/variant_summary.txt.gz" "$SYS"
stamp "restrictions done"

python3 "$C/label_audit.py" "$W/gate/strict/brain_whole/codes" "$W/out/label_audit_whole_pool.json"
cp "$C/genomewide_power.py" "$C/genomewide_rank.py" "$W/gw/codes/"
ln -sf "$ROOT/variant_summary.txt.gz" "$W/gw/variant_summary.txt.gz"
python3 "$W/gw/codes/genomewide_power.py"
python3 "$W/gw/codes/genomewide_rank.py" ARID1A ATM CIC IQSEC2 IRF2BPL KMT2B PPM1D SETD1B SRRM2 ZFHX3
WC="$W/run/brain_whole/codes"
python3 "$C/lcr_counts.py" "$BED" "$WC/gnomad_plof.json" "$WC/derived_v2.json" "$W/out/derived_lcr_counts.json"
python3 "$C/pathogenic_accessions.py" "$ROOT/variant_summary.txt.gz" "$W/pull/raw" "$WC/derived_v2.json" \
  "$WC/gnomad_plof.json" "$W/out/derived_pathogenic_accessions.json"
ARGS=()
for d in "$W"/run/*/; do n=$(basename "$d"); ARGS+=("$n=$W/run/$n"); done
python3 "$C/collate.py" "$W/out/results_summary.json" "${ARGS[@]}"
python3 "$C/assemble_data.py" "$W" "$W/data"
stamp "all done"
