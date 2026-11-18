import gzip
import json
import os
import shutil
import sys

W, D = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
RUN = os.path.join(W, "run")
UNIT_FILES = ["panel_v2.csv", "derived_frame.json", "derived_fsd2.json", "derived_fsd_exclusion.json",
              "derived_fsd_q.json", "derived_labels_submission.json", "derived_labels_submission_strictdual.json",
              "derived_phenotypes.json", "derived_phenotypes_strictdual.json", "derived_pilot.json",
              "derived_pilot_strictdual.json", "derived_pilot_windows.json", "derived_pilot_windows_strictdual.json",
              "derived_pooled_exclusions.json", "derived_pooled_exclusions_strictdual.json", "derived_position_q.json",
              "derived_position_q_strictdual.json", "derived_review.json", "derived_v2.json",
              "resolution_failures.json", "run_panel.json"]
RESTRICTION_FILES = ["derived_frame.json", "derived_fsd2.json", "derived_fsd_q.json", "derived_fsd_exclusion.json",
                     "derived_v2.json", "derived_phenotypes.json", "derived_phenotypes_strictdual.json",
                     "derived_labels_submission.json", "derived_labels_submission_strictdual.json",
                     "derived_position_q.json", "derived_position_q_strictdual.json", "derived_pilot_windows.json",
                     "derived_pilot_windows_strictdual.json", "derived_pooled_exclusions.json",
                     "derived_pooled_exclusions_strictdual.json", "run_panel.json", "resolution_failures.json"]


def put(src, dst, gz=False):
    if not os.path.isfile(src):
        raise SystemExit(f"missing {src}")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if gz:
        with open(src, "rb") as fi, gzip.open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    else:
        shutil.copy(src, dst)


if os.path.exists(D):
    shutil.rmtree(D)
os.makedirs(D)
whole = os.path.join(RUN, "brain_whole", "codes")
for f in UNIT_FILES + ["variant_table.csv"]:
    put(os.path.join(whole, f), os.path.join(D, f))
for f in ("derived_genomewide_power.json", "derived_genomewide_rank.json"):
    put(os.path.join(W, "gw", "codes", f), os.path.join(D, f))
for f in ("derived_lcr_counts.json", "derived_pathogenic_accessions.json", "label_audit_whole_pool.json",
          "gate_brain.json", "grouping_names.json", "results_summary.json"):
    put(os.path.join(W, "out", f), os.path.join(D, f))
put(os.path.join(W, "pull", "union_genes.json"), os.path.join(D, "union_genes.json"))
put(os.path.join(W, "pull", "status.json"), os.path.join(D, "gnomad_pull_status.json"))
put(os.path.join(W, "pull", "no_gene_record.json"), os.path.join(D, "gnomad_no_gene_record.json"))
for u in ("epilepsy", "neurodevelopmental", "neurodegeneration", "malformation"):
    c = os.path.join(RUN, "brain_" + u, "codes")
    for f in UNIT_FILES:
        put(os.path.join(c, f), os.path.join(D, "units", "brain_" + u, f))
    put(os.path.join(c, "variant_table.csv"), os.path.join(D, "units", "brain_" + u, "variant_table.csv.gz"), gz=True)
for s in ("sens_whole_kept", "sens_whole_names10"):
    c = os.path.join(RUN, s, "codes")
    for f in UNIT_FILES:
        put(os.path.join(c, f), os.path.join(D, "sensitivity", s, f))
    put(os.path.join(c, "variant_table.csv"), os.path.join(D, "sensitivity", s, "variant_table.csv.gz"), gz=True)
for run in sorted(os.listdir(RUN)):
    if not run.startswith(("mane_", "lcrsym_", "manelcrsym_", "lcrflag_")) or not os.path.isdir(os.path.join(RUN, run)):
        continue
    c = os.path.join(RUN, run, "codes")
    for f in RESTRICTION_FILES:
        put(os.path.join(c, f), os.path.join(D, "sensitivity", run, f))
    rec = json.load(open(os.path.join(c, "gnomad_plof.json"))).get("restricted")
    json.dump(rec, open(os.path.join(D, "sensitivity", run, "gnomad_restriction.json"), "w"), indent=1)
put(os.path.join(RUN, "vs_lcr_filter.json"), os.path.join(D, "sensitivity", "vs_lcr_filter.json"))
n = sum(len(fn) for _, _, fn in os.walk(D))
print(f"{n} files assembled in {D}")
