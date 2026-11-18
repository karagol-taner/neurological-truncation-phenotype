import argparse, collections, csv, importlib.util, json, os, re

VAGUE = re.compile(r'^(not provided|not specified|see cases|none provided|'
                   r'not applicable|unknown)$', re.I)
FLOOR = 10


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def conditions(gcs):
    cond = collections.defaultdict(set)
    with open(gcs, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = (r.get("AssociatedGenes") or "").strip()
            if s:
                cond[s].add((r.get("DiseaseName") or "").strip())
    return cond


def dual_set(v2, cond, area, excl, rule, sys_inc=None, sys_exc=None):
    out = {}
    iso = {g: e for g, e in v2["genes"].items() if e["tier"] == "isolated"}
    for g, e in iso.items():
        cs = set()
        for s in e["clinvar_symbols"]:
            cs |= cond.get(s, set())
        cs = {c for c in cs if not VAGUE.match(c)}
        outside = sorted(c for c in cs if not (area.search(c) and not excl.search(c)))
        if rule == "strict":
            outside = [c for c in outside if not (sys_inc.search(c) and not sys_exc.search(c))]
            sym = [re.escape(x) for x in e["clinvar_symbols"]]
            umb = re.compile(r'(?<![A-Za-z0-9])(' + "|".join(sym) + r')(?![A-Za-z0-9])', re.I)
            outside = [c for c in outside if not umb.search(c)]
        if outside:
            out[g] = {"n_conditions": len(cs), "other": outside,
                      "positioned_ptv": e["plp_ptv_positioned"]}
    return out, len(iso)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("codes_dir")
    ap.add_argument("--rule", choices=("broad", "strict"), required=True)
    ap.add_argument("--system", help="panel_v2.py whose DISEASE and EXCL define the organ system")
    ap.add_argument("--out", default="derived_phenotypes.json")
    a = ap.parse_args()
    codes = os.path.abspath(a.codes_dir)
    pv = load(os.path.join(codes, "panel_v2.py"), "unit_panel")
    run = os.path.join(codes, "run_panel.json")
    gcs = (json.load(open(run)).get("map_used") if os.path.exists(run) else None) or pv.GCS
    sys_inc = sys_exc = None
    if a.rule == "strict":
        sp = load(a.system or os.path.join(codes, "panel_v2.py"), "system_panel")
        sys_inc, sys_exc = sp.DISEASE, sp.EXCL
    v2 = json.load(open(os.path.join(codes, "derived_v2.json")))
    dual, n_iso = dual_set(v2, conditions(gcs), pv.DISEASE, pv.EXCL, a.rule, sys_inc, sys_exc)
    out = {"dual_phenotype_genes": dual, "n_dual": len(dual),
           "n_dual_ge10": sum(1 for r in dual.values() if r["positioned_ptv"] >= FLOOR),
           "n_isolated": n_iso}
    with open(os.path.join(codes, a.out), "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"rule {a.rule}: isolated {n_iso}, dual {len(dual)}, dual with >= {FLOOR} "
          f"positioned {out['n_dual_ge10']}; map {gcs}")


if __name__ == "__main__":
    main()
