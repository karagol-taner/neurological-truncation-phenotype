import csv, gzip, json, os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
VS = os.path.join(os.path.dirname(HERE), "variant_summary.txt.gz")
PTC = re.compile(r'p\.[A-Za-z]{3}(\d+)(Ter|\*)', re.I)
FS = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs', re.I)

ns = collections.Counter(); fs = collections.Counter()
with gzip.open(VS, "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("Assembly") != "GRCh38":
            continue
        sig = (row.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig:
            continue
        name = row.get("Name") or ""
        if FS.search(name):
            cls = "fs"
        elif PTC.search(name):
            cls = "ns"
        else:
            continue
        for s in (row.get("GeneSymbol") or "").split(";"):
            s = s.strip()
            if not s or s == "-":
                continue
            (fs if cls == "fs" else ns)[s] += 1

genes = sorted(set(ns) | set(fs))
rows = []
for g in genes:
    a, b = fs[g], ns[g]
    t = a + b
    if t:
        rows.append((g, b, a, t, a / t))

def tally(minboth, mintot):
    return [r for r in rows if min(r[1], r[2]) >= minboth and r[3] >= mintot]

print(f"genes with at least one pathogenic truncating variant: {len(rows)}")
for mb, mt in ((1, 10), (3, 10), (5, 20), (5, 30), (10, 40)):
    s = tally(mb, mt)
    print(f"  both classes >= {mb} and total >= {mt}: {len(s)} genes")

pool = [r for r in rows if r[3] >= 20]
print(f"\ngenes with >= 20 pathogenic truncating variants: {len(pool)}")
shares = sorted(r[4] for r in pool)
q = lambda p: round(shares[int(p * (len(shares) - 1))], 3)
print(f"\nframeshift share across the pool: min {q(0)}, 10th {q(.1)}, median {q(.5)}, "
      f"90th {q(.9)}, max {q(1)}")
json.dump({"n_genes_with_ptv": len(rows),
           "qualifying": {f"both_ge{mb}_total_ge{mt}": len(tally(mb, mt))
                          for mb, mt in ((1,10),(3,10),(5,20),(5,30),(10,40))},
           "pool_ge20": len(pool),
           "share_quantiles": {"min": q(0), "p10": q(.1), "median": q(.5), "p90": q(.9), "max": q(1)}},
          open(os.path.join(HERE, "derived_genomewide_power.json"), "w"), indent=1, sort_keys=True)
