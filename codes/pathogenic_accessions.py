import collections, csv, glob, gzip, json, os, re, sys

vs, raw, v2p, gnp, out = sys.argv[1:6]
mane = {}
for f in glob.glob(os.path.join(raw, "gene_*.json.gz")):
    g = ((json.load(gzip.open(f, "rt")).get("response") or {}).get("data") or {}).get("gene")
    if g and g.get("mane_select_transcript"):
        mane[g["symbol"]] = g["mane_select_transcript"]["refseq_id"]
v2 = json.load(open(v2p))["genes"]
iso = {g for g, e in v2.items() if e["tier"] == "isolated"}
sym2gene = {s: g for g, e in v2.items() for s in e["clinvar_symbols"]}
opener = gzip.open if gnp.endswith(".gz") else open
gn = json.load(opener(gnp, "rt"))["genes"]
gnomad_symbol = {(r.get("clinvar_symbol") or s): s for s, r in gn.items()}
PTC = re.compile(r"p\.[A-Za-z]{3}(\d+)(Ter|\*)", re.I)
FS = re.compile(r"p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs", re.I)
NAME = re.compile(r"^(N[MR]_\d+)\.\d+\(([^)]+)\)")
total, per = collections.Counter(), collections.defaultdict(collections.Counter)
with gzip.open(vs, "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("Assembly") != "GRCh38":
            continue
        hits = {sym2gene[s] for s in (row.get("GeneSymbol") or "").split(";")
                if s in sym2gene and sym2gene[s] in iso}
        if not hits:
            continue
        sig = (row.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig:
            continue
        name = row.get("Name") or ""
        if not (FS.search(name) or PTC.search(name)):
            continue
        g = sorted(hits)[0]
        m = NAME.match(name)
        ref = mane.get(gnomad_symbol.get(g, g))
        k = ("no accession" if not m else "no MANE Select record" if ref is None
             else "MANE Select" if m.group(1) == ref else "other transcript")
        total[k] += 1
        per[g][k] += 1
json.dump({"counts": dict(total), "per_gene": {g: dict(c) for g, c in sorted(per.items())},
           "mane_select_refseq": dict(sorted(mane.items()))}, open(out, "w"), indent=1)
print(dict(total), sum(total.values()))
