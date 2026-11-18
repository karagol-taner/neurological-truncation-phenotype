import argparse, collections, csv, json, os, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import dual_phenotypes as dp

UNITS = {"brain_whole": os.path.join(HERE, "panel_v2.py")}
for u in ("brain_epilepsy", "brain_neurodevelopmental", "brain_neurodegeneration", "brain_malformation"):
    UNITS[u] = os.path.join(ROOT, "units", u, "codes", "panel_v2.py")
CLINVAR = ("variant_summary.txt.gz", "submission_summary.txt.gz", "gene_condition_source_id")
AMEND_FROM = "|osteo|chondro', re.I)"
AMEND_TO = "|osteo|chondro|lung cysts|wilms tumor', re.I)"
VARIANTS = {
    "screen": ("first", 0, False, "names"),
    "first_dicer": ("first", 0, True, "names"),
    "strict": ("strict", 0, False, "names"),
    "strict_dicer": ("strict", 0, True, "names"),
    "strict_dicer_g20": ("strict", 20, True, "names"),
    "strict_dicer_g10": ("strict", 10, True, "names"),
    "strict_dicer_r20": ("strict", 20, True, "routes"),
    "strict_dicer_r10": ("strict", 10, True, "routes"),
}
FLOOR = 10


def amended(text):
    line = [l for l in text.splitlines() if l.startswith("SYNDROMIC = ")]
    assert len(line) == 1 and line[0].endswith(AMEND_FROM), "SYNDROMIC line not as expected"
    return text.replace(line[0], line[0][:-len(AMEND_FROM)] + AMEND_TO)


def prepare(work, variant, unit):
    res, grp, dicer, scope = VARIANTS[variant]
    d = os.path.join(work, variant, unit)
    os.makedirs(os.path.join(d, "codes"), exist_ok=True)
    src = open(UNITS[unit], encoding="utf-8").read().replace(AMEND_TO, AMEND_FROM)
    open(os.path.join(d, "codes", "panel_v2.py"), "w", encoding="utf-8").write(
        amended(src) if dicer else src)
    for f in CLINVAR:
        link = os.path.join(d, f)
        if not os.path.lexists(link):
            os.symlink(os.path.realpath(os.path.join(ROOT, f)), link)
    return d


def run_one(args):
    work, variant, unit = args
    res, grp, dicer, scope = VARIANTS[variant]
    d = prepare(work, variant, unit)
    codes = os.path.join(d, "codes")
    if not os.path.exists(os.path.join(codes, "derived_v2.json")):
        with open(os.path.join(d, "run.log"), "w") as log:
            r = subprocess.run([sys.executable, os.path.join(HERE, "run_panel.py"), codes,
                                "--resolution", res, "--grouping", str(grp),
                                "--grouping-scope", scope],
                               stdout=log, stderr=subprocess.STDOUT)
        if r.returncode:
            raise SystemExit(f"run_panel failed for {variant}/{unit}, see {d}/run.log")
    for rule, out in (("broad", "derived_phenotypes.json"), ("strict", "derived_phenotypes_strict.json")):
        r = subprocess.run([sys.executable, os.path.join(HERE, "dual_phenotypes.py"), codes,
                            "--rule", rule, "--system", UNITS["brain_whole"], "--out", out],
                           capture_output=True, text=True)
        if r.returncode:
            raise SystemExit(f"dual_phenotypes failed for {variant}/{unit}: {r.stderr}")
    return variant, unit


def summarise(work, variant, unit, full_cond):
    codes = os.path.join(work, variant, unit, "codes")
    v2 = json.load(open(os.path.join(codes, "derived_v2.json")))
    run = json.load(open(os.path.join(codes, "run_panel.json")))
    ph = json.load(open(os.path.join(codes, "derived_phenotypes.json")))
    phs = json.load(open(os.path.join(codes, "derived_phenotypes_strict.json")))
    genes = v2["genes"]
    iso = {g: e for g, e in genes.items() if e["tier"] == "isolated"}
    frame = {g for g, e in iso.items() if e["plp_ptv_positioned"] >= FLOOR}
    both = {g for g in frame if iso[g]["plp_nonsense"] >= 1 and iso[g]["plp_frameshift"] >= 1}
    pos_d = {g for g, r in ph["dual_phenotype_genes"].items() if r["positioned_ptv"] >= FLOOR}
    pos_s = {g for g, r in phs["dual_phenotype_genes"].items() if r["positioned_ptv"] >= FLOOR}
    ranked = sorted(iso, key=lambda g: (-iso[g]["plp_ptv"], g))
    top = ranked[0] if ranked else None
    tot = v2["isolated_ptv_total"]
    rec = {
        "resolution": run["resolution"], "grouping_threshold": run["grouping_threshold"],
        "grouping_scope": run.get("grouping_scope", "names" if run["grouping_threshold"] else None),
        "grouping_names": run.get("grouping_names", run.get("grouping_names_dropped", 0)),
        "symbols_whose_only_route_is_a_grouping_name": run.get("symbols_whose_only_route_is_a_grouping_name"),
        "rows_grch38": v2["clinvar_rows_matched_grch38"],
        "symbols": v2["panel_symbols_matched"], "genes": v2["panel_genes_after_alias_merge"],
        "isolated": v2["n_isolated"], "syndromic_only": v2["n_syndromic_only"],
        "alias_merges": v2["alias_merges"], "resolution_failures": run["resolution_failures"],
        "isolated_ptv_total": tot,
        "frame_pool": len(frame), "frame_both_classes": len(both),
        "position_pool_broad_rule": len(pos_d), "position_pool_strict_rule": len(pos_s),
        "gate4_both_and_dual_broad_rule": len(both & pos_d),
        "gate4_both_and_dual_strict_rule": len(both & pos_s),
        "largest_gene": top, "largest_gene_ptv": iso[top]["plp_ptv"] if top else None,
        "largest_gene_share": round(iso[top]["plp_ptv"] / tot, 3) if top and tot else None,
        "five_largest": ranked[:5],
        "without_ATM": {"frame_pool": len(frame - {"ATM"}), "frame_both_classes": len(both - {"ATM"}),
                        "position_pool_broad_rule": len(pos_d - {"ATM"}),
                        "position_pool_strict_rule": len(pos_s - {"ATM"}),
                        "ATM_in_isolated_tier": "ATM" in iso},
        "DICER1_tier": genes.get("DICER1", {}).get("tier"),
        "isolated_genes": sorted(iso), "frame_genes": sorted(frame),
        "position_genes_broad_rule": sorted(pos_d), "position_genes_strict_rule": sorted(pos_s),
    }
    if run["grouping_threshold"]:

        sysm = dp.load(UNITS["brain_whole"], "sys_whole")
        um = dp.load(os.path.join(codes, "panel_v2.py"), "unit_" + variant + unit)
        dd, _ = dp.dual_set(v2, full_cond, um.DISEASE, um.EXCL, "broad")
        ds, _ = dp.dual_set(v2, full_cond, um.DISEASE, um.EXCL, "strict", sysm.DISEASE, sysm.EXCL)
        rec["position_pool_broad_rule_second_conditions_unfiltered"] = sum(
            1 for r in dd.values() if r["positioned_ptv"] >= FLOOR)
        rec["position_pool_strict_rule_second_conditions_unfiltered"] = sum(
            1 for r in ds.values() if r["positioned_ptv"] >= FLOOR)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--variants", nargs="*", default=list(VARIANTS))
    a = ap.parse_args()
    work = os.path.abspath(a.work)
    jobs = [(work, v, u) for v in a.variants for u in UNITS]
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for v, u in ex.map(run_one, jobs):
            print("done", v, u, flush=True)
    full_cond = dp.conditions(os.path.join(ROOT, "gene_condition_source_id"))
    out = {"floor": FLOOR, "variants": {v: {"resolution": VARIANTS[v][0], "grouping_threshold": VARIANTS[v][1],
                                            "dicer1_amendment": VARIANTS[v][2], "grouping_scope": VARIANTS[v][3]}
                                        for v in VARIANTS},
           "units": {}}
    prev = os.path.join(HERE, "gate_brain.json")
    if os.path.exists(prev):
        out["units"] = json.load(open(prev)).get("units", {})
    for v in a.variants:
        for u in UNITS:
            out["units"].setdefault(v, {})[u] = summarise(work, v, u, full_cond)
    json.dump(out, open(prev, "w"), indent=1, sort_keys=True)
    for v in a.variants:
        for u in UNITS:
            r = out["units"][v][u]
            print(f"{v:18s} {u:26s} genes {r['genes']:5d} iso {r['isolated']:4d} syn {r['syndromic_only']:4d} "
                  f"frame {r['frame_pool']:4d} both {r['frame_both_classes']:4d} pos_d {r['position_pool_broad_rule']:4d} "
                  f"pos_s {r['position_pool_strict_rule']:4d} top {r['largest_gene']} {r['largest_gene_share']} "
                  f"DICER1 {r['DICER1_tier']}")


if __name__ == "__main__":
    main()
