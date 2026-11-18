import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "genomewide_power.py")).read()
ns = {"__file__": os.path.join(HERE, "genomewide_power.py")}
exec(compile(src[:src.index("def tally(")], "genomewide_power.table", "exec"), ns)
rows = {r[0]: r for r in ns["rows"]}
pool = [r for r in ns["rows"] if r[3] >= 20]
out = {"pool_ge20": len(pool), "n_genes_with_ptv": len(ns["rows"]), "genes": {}}
for g in sys.argv[1:]:
    r = rows.get(g)
    if not r:
        out["genes"][g] = None
        continue
    out["genes"][g] = {"nonsense": r[1], "frameshift": r[2], "total": r[3], "share": round(r[4], 3),
                       "in_pool_ge20": r[3] >= 20,
                       "n_pool_at_or_above_share": sum(1 for x in pool if x[4] >= r[4]),
                       "n_pool_at_or_below_share": sum(1 for x in pool if x[4] <= r[4])}
json.dump(out, open(os.path.join(HERE, "derived_genomewide_rank.json"), "w"), indent=1, sort_keys=True)
for g, r in out["genes"].items():
    print(g, r)
