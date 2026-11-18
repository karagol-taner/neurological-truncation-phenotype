import csv, json, os, statistics
from scipy.stats import mannwhitneyu, fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
FIXED = "ATM"
THRESH = (10, 20, 40, 80)


def median(xs):
    if not xs:
        return None
    m = statistics.median(xs)
    return int(m) if float(m).is_integer() else m


def compare(p, q):
    rec = {"n_path": len(p), "n_pop": len(q), "med_path": median(p), "med_pop": median(q)}
    if len(p) >= 3 and len(q) >= 3:
        rec["mwu_p"] = float(mannwhitneyu(p, q, alternative="two-sided")[1])
    for t in THRESH:
        a = sum(1 for x in p if x >= t)
        c = sum(1 for x in q if x >= t)
        rec["ge%d" % t] = {"path": "%d/%d" % (a, len(p)), "pop": "%d/%d" % (c, len(q)),
                           "p": float(fisher_exact([[a, len(p) - a], [c, len(q) - c]])[1])}
    return rec


def neo_tail():
    path = os.path.join(HERE, "variant_table.csv")
    if not os.path.exists(path):
        return None
    vt = list(csv.DictReader(open(path)))
    rows = [(r["gene"], r["source"], int(r["fsd"])) for r in vt
            if r["variant_class"] == "frameshift" and r["fsd"] not in ("", None)]

    def tails(source, drop=None):
        return [x for g, s, x in rows if s == source and g != drop]

    out = {"pooled_all": compare(tails("clinvar_pathogenic"), tails("gnomad_plof")),
           f"pooled_excluding_{FIXED}": compare(tails("clinvar_pathogenic", FIXED),
                                                tails("gnomad_plof", FIXED))}
    genes = sorted({g for g, s, x in rows if s == "clinvar_pathogenic"})
    loo = {}
    for g in genes:
        p, q = tails("clinvar_pathogenic", g), tails("gnomad_plof", g)
        if len(p) >= 3 and len(q) >= 3:
            loo[g] = float(mannwhitneyu(p, q, alternative="two-sided")[1])
    if loo:
        carrier = max(loo, key=lambda g: (loo[g], g))
        out["carrying_gene"] = carrier
        out["carrying_gene_rule"] = "the gene whose removal raises the pooled Mann-Whitney P value most"
        out[f"pooled_excluding_{carrier}"] = compare(tails("clinvar_pathogenic", carrier),
                                                     tails("gnomad_plof", carrier))
        out["leave_one_out_mwu_p"] = loo
    return out


def windows(fname):
    path = os.path.join(HERE, fname)
    if not os.path.exists(path):
        return None
    d = json.load(open(path))
    res = {}
    defs = sorted({w for r in d["genes"].values() for w in r["windows"]})
    for w in defs:
        tabs = {g: r["windows"][w] for g, r in d["genes"].items() if w in r["windows"]
                and (r["windows"][w]["area_in"] + r["windows"][w]["area_out"])
                and (r["windows"][w]["other_in"] + r["windows"][w]["other_out"])}

        def pooled(drop=None):
            a = sum(t["area_in"] for g, t in tabs.items() if g != drop)
            b = sum(t["area_out"] for g, t in tabs.items() if g != drop)
            c = sum(t["other_in"] for g, t in tabs.items() if g != drop)
            dd = sum(t["other_out"] for g, t in tabs.items() if g != drop)
            return {"area_in": a, "area_out": b, "other_in": c, "other_out": dd,
                    "n_genes": len([g for g in tabs if g != drop]),
                    "p": float(fisher_exact([[a, b], [c, dd]])[1])}

        rec = {"pooled_all": pooled(), f"pooled_excluding_{FIXED}": pooled(FIXED)}
        loo = {g: pooled(g)["p"] for g in tabs}
        if loo:
            carrier = max(loo, key=lambda g: (loo[g], g))
            rec["carrying_gene"] = carrier
            rec[f"pooled_excluding_{carrier}"] = pooled(carrier)
            rec["leave_one_out_p"] = loo
        res[w] = rec
    return res


out = {"fixed_exclusion": FIXED,
       "neo_tail": neo_tail(),
       "windows_submission_level": windows("derived_labels_submission.json"),
       "windows_aggregated_list": windows("derived_pilot_windows.json")}
json.dump(out, open(os.path.join(HERE, "derived_pooled_exclusions.json"), "w"), indent=1, sort_keys=True)
nt = out["neo_tail"]
if nt:
    for k in sorted(k for k in nt if k.startswith("pooled")):
        r = nt[k]
        print(f"neo-tail {k:28s} n {r['n_path']}/{r['n_pop']} medians {r['med_path']}/{r['med_pop']} "
              f"p {r.get('mwu_p')}")
for src in ("windows_submission_level", "windows_aggregated_list"):
    if out[src]:
        for w, rec in out[src].items():
            for k in sorted(k for k in rec if k.startswith("pooled")):
                print(f"{src} {w:14s} {k:28s} p {rec[k]['p']:.3g} genes {rec[k]['n_genes']}")
