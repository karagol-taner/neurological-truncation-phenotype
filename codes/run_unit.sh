#!/bin/bash
set -euo pipefail
U=$(cd "$1" && pwd); G=$2; SCOPE=$3; SYS=$4; C=$U/codes; L=$U/logs
mkdir -p "$L"; cd "$C"
step() { local n=$1; shift; echo "$(date -u +%H:%M:%S) $n" >> "$L/steps.txt"; "$@" > "$L/$n.txt" 2>&1 || { echo "FAILED $n" >> "$L/steps.txt"; exit 2; }; }
[ -s gnomad_plof.json ] || { echo "FAILED gnomad_plof.json missing" >> "$L/steps.txt"; exit 2; }
step 01_panel python3 run_panel.py "$C" --resolution strict --grouping "$G" --grouping-scope "$SCOPE"
step 02_review_filter python3 review_filter.py
step 03_fsd2 python3 fsd2.py
step 04_export_variant_table python3 export_variant_table.py
step 05_frame_dependence python3 frame_dependence.py
step 06_fsd_summary python3 fsd_summary.py
step 06b_fsd_q python3 fsd_q.py
POS="derived_pilot.json derived_pilot_windows.json derived_labels_submission.json derived_position_q.json derived_pooled_exclusions.json"
step 07_dual_strict python3 dual_phenotypes.py "$C" --rule strict --system "$SYS"
step 08_pilot_positions_strictdual python3 pilot_positions.py
step 09_pilot_windows_strictdual python3 pilot_windows.py
step 10_labels_submission_level_strictdual python3 labels_submission_level.py
step 11_position_q_strictdual python3 position_q.py
step 12_pooled_exclusions_strictdual python3 pooled_exclusions.py
mv derived_phenotypes.json derived_phenotypes_strictdual.json
for f in $POS; do mv "$f" "${f%.json}_strictdual.json"; done
step 13_dual_broad python3 dual_phenotypes.py "$C" --rule broad
step 14_pilot_positions python3 pilot_positions.py
step 15_pilot_windows python3 pilot_windows.py
step 16_labels_submission_level python3 labels_submission_level.py
step 17_position_q python3 position_q.py
step 18_pooled_exclusions python3 pooled_exclusions.py
echo "$(date -u +%H:%M:%S) DONE" >> "$L/steps.txt"
