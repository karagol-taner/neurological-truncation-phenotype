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


d = json.load(open(os.path.join(HERE, "derived_fsd2.json")))
out = {}
pairs = [(g, r["mwu_p"]) for g, r in d["genes"].items() if r.get("mwu_p") is not None]
q = bh(pairs)
out["mwu"] = {"n_genes_tested": len(pairs), "genes": {g: {"p": p, "q": q[g]} for g, p in sorted(pairs)},
              "q_lt_0.05": sorted(g for g in q if q[g] < 0.05)}
for t in d["thresholds"]:
    k = f"ge{t}"
    pairs = [(g, r["long"][k]["p"]) for g, r in d["genes"].items() if k in r.get("long", {})]
    q = bh(pairs)
    out[f"long_{k}"] = {"n_genes_tested": len(pairs),
                        "genes": {g: {"p": p, "q": q[g]} for g, p in sorted(pairs)},
                        "q_lt_0.05": sorted(g for g in q if q[g] < 0.05)}
json.dump(out, open(os.path.join(HERE, "derived_fsd_q.json"), "w"), indent=1, sort_keys=True)
for k, r in out.items():
    print(f"{k:10s} tested {r['n_genes_tested']:4d} q<0.05 {len(r['q_lt_0.05'])}: {r['q_lt_0.05'][:30]}")
