#!/bin/bash
set -euo pipefail
R=$(cd "$1" && pwd); BED=$2; VS=$3; SYS=$4
C=$(cd "$(dirname "$0")" && pwd)
python3 "$C/vs_lcr_filter.py" "$BED" "$VS" "$R/variant_summary_outside_lcr.txt.gz" "$R/vs_lcr_filter.json"
G=$R/brain_whole/codes/gnomad_plof.json
python3 "$C/gnomad_restrict.py" "$G" "$R/gnomad_mane.json" --mane
python3 "$C/gnomad_restrict.py" "$G" "$R/gnomad_lcr.json" --lcr-bed "$BED"
python3 "$C/gnomad_restrict.py" "$G" "$R/gnomad_mane_lcr.json" --mane --lcr-bed "$BED"
python3 "$C/gnomad_restrict.py" "$G" "$R/gnomad_flag.json" --lcr-flag
setup() {  # setup NEW_RUN FROM_UNIT GNOMAD_FILE OUTSIDE_LCR(0/1)
  rm -rf "$R/$1"; cp -a "$R/$2" "$R/$1"; rm -rf "$R/$1/logs"
  ( cd "$R/$1/codes" && rm -f derived_*.json panel_v2.csv resolution_failures.json run_panel.json \
      variant_table.csv gene_condition_source_id.grouping_* )
  cp "$3" "$R/$1/codes/gnomad_plof.json"
  if [ "$4" = 1 ]; then rm -f "$R/$1/variant_summary.txt.gz"; ln -s "$R/variant_summary_outside_lcr.txt.gz" "$R/$1/variant_summary.txt.gz"; fi
}
for U in whole epilepsy neurodevelopmental neurodegeneration malformation; do
  setup mane_$U brain_$U "$R/gnomad_mane.json" 0
  setup lcrsym_$U brain_$U "$R/gnomad_lcr.json" 1
  setup manelcrsym_$U brain_$U "$R/gnomad_mane_lcr.json" 1
done
setup lcrflag_whole brain_whole "$R/gnomad_flag.json" 0
for D in "$R"/mane_* "$R"/lcrsym_* "$R"/manelcrsym_* "$R"/lcrflag_whole; do
  bash "$D/codes/run_unit.sh" "$D" 10 routes "$SYS"
done
