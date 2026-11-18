#!/usr/bin/env python3

import csv, json, os, statistics, sys
from scipy.stats import mannwhitneyu, fisher_exact

def median(xs):

    if not xs:
        return None
    m = statistics.median(xs)
    return int(m) if float(m).is_integer() else m

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = HERE if os.path.exists(os.path.join(HERE, "variant_table.csv")) \
       else os.path.join(HERE, "..", "data")
THRESH = (10, 20, 40, 80)

with open(os.path.join(DATA, "variant_table.csv")) as fh:
    VT = list(csv.DictReader(fh))

def tails(source, gene=None, drop=None):
    return [int(r["fsd"]) for r in VT
            if r["source"] == source and r["variant_class"] == "frameshift"
            and r["fsd"] not in ("", None)
            and (gene is None or r["gene"] == gene)
            and (drop is None or r["gene"] != drop)]

def compare(p, q):
    rec = {"n_path": len(p), "n_pop": len(q),
           "med_path": median(p),
           "med_pop": median(q)}
    if len(p) >= 3 and len(q) >= 3:
        rec["mwu_p"] = float(mannwhitneyu(p, q, alternative="two-sided")[1])
    for t in THRESH:
        a = sum(1 for x in p if x >= t)
        c = sum(1 for x in q if x >= t)
        rec["ge%d" % t] = {"path": "%d/%d" % (a, len(p)), "pop": "%d/%d" % (c, len(q)),
                           "p": float(fisher_exact([[a, len(p) - a], [c, len(q) - c]])[1])}
    return rec

ex = {"pooled_all": compare(tails("clinvar_pathogenic"), tails("gnomad_plof")),
      "pooled_excluding_ATM": compare(tails("clinvar_pathogenic", drop="ATM"),
                                      tails("gnomad_plof", drop="ATM"))}
json.dump(ex, open(os.path.join(DATA, "derived_fsd_exclusion.json"), "w"),
          indent=1, sort_keys=True)

path = os.path.join(DATA, "derived_fsd2.json")
fs2 = json.load(open(path))
changed = []
for g, rec in fs2["genes"].items():
    for key, source in (("median_fsd_pathogenic", "clinvar_pathogenic"),
                        ("median_fsd_population", "gnomad_plof")):
        xs = tails(source, gene=g)
        new, old = median(xs), rec.get(key)
        if repr(old) != repr(new):
            changed.append((g, key, old, new))
        rec[key] = new
allp, allq = tails("clinvar_pathogenic"), tails("gnomad_plof")
for key, xs in (("median_pathogenic", allp), ("median_population", allq)):
    new, old = median(xs), fs2.get("pooled", {}).get(key)
    if repr(old) != repr(new):
        changed.append(("pooled", key, old, new))
    fs2.setdefault("pooled", {})[key] = new
json.dump(fs2, open(path, "w"), indent=1, sort_keys=True)

print("wrote derived_fsd_exclusion.json and refreshed derived_fsd2.json")
for g, key, old, new in changed:
    print("   %-10s %-24s %s -> %s" % (g, key, old, new))
if not changed:
    print("   medians already agreed")
