import json, os

HERE = os.path.dirname(os.path.abspath(__file__))


def bh(pairs):

    xs = sorted(pairs, key=lambda kv: kv[1])
    n = len(xs)
    q, prev = {}, 1.0
    for i in range(n - 1, -1, -1):
        k, p = xs[i]
        prev = min(prev, p * n / (i + 1))
        q[k] = round(min(prev, 1.0), 6)
    return q


def windows_q(path):
    d = json.load(open(path))
    out = {}
    defs = sorted({w for r in d["genes"].values() for w in r["windows"]})
    for w in defs:
        for comp in ("p_area_vs_other", "p_area_vs_gnomad"):
            pairs = [(g, r["windows"][w][comp]) for g, r in d["genes"].items()
                     if w in r["windows"] and r["windows"][w].get(comp) is not None]
            q = bh(pairs)
            out.setdefault(w, {})[comp] = {
                "n_genes_tested": len(pairs),
                "genes": {g: {"p": p, "q": q[g]} for g, p in sorted(pairs)},
                "q_lt_0.05": sorted(g for g in q if q[g] < 0.05)}
    return out


res = {}
for name, f in (("submission_level", "derived_labels_submission.json"),
                ("aggregated_list", "derived_pilot_windows.json")):
    p = os.path.join(HERE, f)
    if os.path.exists(p):
        res[name] = windows_q(p)
p = os.path.join(HERE, "derived_pilot.json")
if os.path.exists(p):
    d = json.load(open(p))
    for comp in ("mwu_p_area_vs_other", "mwu_p_area_vs_gnomad"):
        pairs = [(g, r[comp]) for g, r in d["genes"].items() if r.get(comp) is not None]
        q = bh(pairs)
        res.setdefault("median_test", {})[comp] = {
            "n_genes_tested": len(pairs),
            "genes": {g: {"p": pv, "q": q[g]} for g, pv in sorted(pairs)},
            "q_lt_0.05": sorted(g for g in q if q[g] < 0.05)}
json.dump(res, open(os.path.join(HERE, "derived_position_q.json"), "w"), indent=1, sort_keys=True)
for src, defs in res.items():
    for w, comps in defs.items():
        if "n_genes_tested" in comps:
            print(f"{src:16s} {w:28s} tested {comps['n_genes_tested']:4d} q<0.05 {comps['q_lt_0.05']}")
            continue
        for comp, r in comps.items():
            print(f"{src:16s} {w:14s} {comp:24s} tested {r['n_genes_tested']:4d} q<0.05 {r['q_lt_0.05']}")
