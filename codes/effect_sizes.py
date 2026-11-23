import csv
import json
import math
import os
import sys

from scipy.stats import fisher_exact, mannwhitneyu, norm

D = sys.argv[1]
Z = float(norm.ppf(0.975))


def load(name):
    return json.load(open(os.path.join(D, name)))


def die(msg):
    sys.exit("FATAL effect_sizes.py: " + msg)


def same(a, b):
    return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-300)


def wilson(x, n):
    p = x / n
    den = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / den
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / den
    return c - h, c + h


def newcombe(x1, n1, x2, n2):
    p1, p2 = x1 / n1, x2 / n2
    l1, u1 = wilson(x1, n1)
    l2, u2 = wilson(x2, n2)
    d = p1 - p2
    return d, d - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), d + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)


def hodges_lehmann(x, y):
    ds = sorted(a - b for a in x for b in y)
    n, m = len(x), len(y)
    N = n * m
    k = max(1, int(math.floor(N / 2 - Z * math.sqrt(n * m * (n + m + 1) / 12))))
    return (ds[(N - 1) // 2] + ds[N // 2]) / 2, ds[k - 1], ds[N - k]


W = load("results_summary.json")["brain_whole"]
out = {"confidence": 0.95, "window_nmd_50nt_rule": {}, "frame_composition": {},
       "frame_recurrent": {}, "tail_mwu": {}}

LS = load("derived_labels_submission.json")["genes"]
PQ = load("derived_position_q.json")["submission_level"]["nmd_50nt_rule"]["p_area_vs_gnomad"]["genes"]
for g in W["position_broad_rule"]["windows"]["nmd_50nt_rule"]["p_area_vs_gnomad"]["q_lt_0.05"]:
    w = LS[g]["windows"]["nmd_50nt_rule"]
    a, b, e, f = w["area_in"], w["area_out"], w["gnomad_in"], w["gnomad_out"]
    if a + b != LS[g]["n_area"] or e + f != LS[g]["n_gnomad_hc_pass"]:
        die(g + ": window counts do not add up to the gene totals")
    if not same(fisher_exact([[a, b], [e, f]])[1], PQ[g]["p"]):
        die(g + ": the window counts do not reproduce the tested P value")
    d, lo, hi = newcombe(a, a + b, e, e + f)
    out["window_nmd_50nt_rule"][g] = {"pathogenic_in": a, "pathogenic_n": a + b, "population_in": e,
                                     "population_n": e + f, "difference": d, "ci": [lo, hi]}

FR = load("derived_frame.json")["genes"]
sig = W["frame"]["q_lt_0.05"]
for g in sorted(set(sig["T1"]) | set(sig["T4"])):
    r = FR[g]
    pns, pfs, gns, gfs = r["n_path_nonsense"], r["n_path_frameshift"], r["n_gnomad_nonsense"], r["n_gnomad_frameshift"]
    rns, rfs = r["n_gnomad_recurrent_nonsense"], r["n_gnomad_recurrent_frameshift"]
    if not same(fisher_exact([[pns, pfs], [gns, gfs]])[1], r["T1_p"]):
        die(g + ": the frame counts do not reproduce the composition P value")
    if not same(fisher_exact([[pns, pfs], [rns, rfs]])[1], r["T4_p"]):
        die(g + ": the frame counts do not reproduce the recurrent P value")
    for key, x2, n2 in (("frame_composition", gfs, gfs + gns), ("frame_recurrent", rfs, rfs + rns)):
        d, lo, hi = newcombe(pfs, pfs + pns, x2, n2)
        out[key][g] = {"pathogenic_frameshift": pfs, "pathogenic_n": pfs + pns, "population_frameshift": x2,
                       "population_n": n2, "difference": d, "ci": [lo, hi]}

F2 = load("derived_fsd2.json")["genes"]
tail_genes = load("derived_fsd_q.json")["mwu"]["q_lt_0.05"]
if sorted(tail_genes) != sorted(W["frame"]["neo_tail_q_lt_0.05"]["mwu"]):
    die("the tail genes of derived_fsd_q.json and results_summary.json differ")
dist = {g: {"clinvar_pathogenic": [], "gnomad_plof": []} for g in tail_genes}
with open(os.path.join(D, "variant_table.csv"), newline="") as fh:
    for row in csv.DictReader(fh):
        if row["gene"] in dist and row["variant_class"] == "frameshift" and row["fsd"] != "":
            dist[row["gene"]][row["source"]].append(float(row["fsd"]))
for g in tail_genes:
    x, y = dist[g]["clinvar_pathogenic"], dist[g]["gnomad_plof"]
    if len(x) != F2[g]["n_pathogenic"] or len(y) != F2[g]["n_population"] or F2[g]["n_pathogenic_censored"]:
        die(g + ": the per-variant table does not hold the tested tails")
    if not same(mannwhitneyu(x, y, alternative="two-sided").pvalue, F2[g]["mwu_p"]):
        die(g + ": the per-variant tails do not reproduce the tested P value")
    s, lo, hi = hodges_lehmann(x, y)
    out["tail_mwu"][g] = {"n_pathogenic": len(x), "n_population": len(y), "shift": s, "ci": [lo, hi]}

json.dump(out, open(os.path.join(D, "derived_effect_sizes.json"), "w"), indent=1, sort_keys=True)
print("wrote derived_effect_sizes.json: %d window, %d frame, %d tail genes"
      % (len(out["window_nmd_50nt_rule"]), len(out["frame_composition"]), len(out["tail_mwu"])))
