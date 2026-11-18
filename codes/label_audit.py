import csv, gzip, json, os, re, sys, collections

codes = os.path.abspath(sys.argv[1])
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, codes)
from area import IN_AREA
SCRATCH = os.path.dirname(codes)
VS = os.path.join(SCRATCH, "variant_summary.txt.gz")
SS = os.path.join(SCRATCH, "submission_summary.txt.gz")
PH = json.load(open(os.path.join(codes, "derived_phenotypes.json")))
V2 = json.load(open(os.path.join(codes, "derived_v2.json")))
VAGUE = re.compile(r'not provided|not specified|see cases|none provided|not applicable|unknown|'
                   r'^-$|related disorder', re.I)
PTC = re.compile(r'p\.[A-Za-z]{3}(\d+)(Ter|\*)', re.I)
FS = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs', re.I)
POOL = sorted(g for g, r in PH["dual_phenotype_genes"].items() if r["positioned_ptv"] >= 10)
sym2gene = {s: g for g, e in V2["genes"].items() for s in e["clinvar_symbols"]}

vid = {}
with gzip.open(VS, "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("Assembly") != "GRCh38": continue
        hits = {sym2gene[s] for s in (row.get("GeneSymbol") or "").split(";")
                if s in sym2gene and sym2gene[s] in POOL}
        if not hits: continue
        sig = (row.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig: continue
        m = PTC.search(row.get("Name") or "") or FS.search(row.get("Name") or "")
        if not m: continue
        vid[row["VariationID"]] = sorted(hits)[0]

side_terms = {"area": collections.Counter(), "other": collections.Counter(), "vague": collections.Counter()}
per_var = collections.defaultdict(collections.Counter)
gene_terms = collections.defaultdict(collections.Counter)
with gzip.open(SS, "rt", encoding="utf-8", errors="replace") as fh:
    header = None
    for line in fh:
        if line.startswith("#"):
            if "VariationID" in line:
                header = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if header is None: continue
        f = line.rstrip("\n").split("\t")
        if len(f) < len(header): continue
        r = dict(zip(header, f))
        v = r.get("VariationID")
        if v not in vid: continue
        sig = (r.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig: continue
        ph = re.sub(r'^[A-Z0-9:]+[:;]', '', r.get("ReportedPhenotypeInfo") or "")
        if VAGUE.search(ph) and not IN_AREA.search(ph):
            b = "vague"
        else:
            b = "area" if IN_AREA.search(ph) else "other"
        side_terms[b][ph] += 1
        per_var[v][b] += 1
        gene_terms[vid[v]][(b, ph)] += 1

buckets = collections.Counter()
for v, g in vid.items():
    c = per_var[v]
    if c["area"] and not c["other"]: buckets["area"] += 1
    elif c["other"] and not c["area"]: buckets["other"] += 1
    elif c["area"] and c["other"]: buckets["both"] += 1
    else: buckets["unlabelled"] += 1
out = {"pool": POOL, "n_variants": len(vid), "variant_buckets": dict(buckets),
       "submissions_by_side": {k: sum(v.values()) for k, v in side_terms.items()},
       "top_terms": {k: v.most_common(60) for k, v in side_terms.items()},
       "per_gene_top_other": {g: [(t, n) for (b, t), n in c.most_common() if b == "other"][:8]
                                   for g, c in gene_terms.items()}}
json.dump(out, open(sys.argv[2], "w"), indent=1)
print("pool", len(POOL), "variants", len(vid), "buckets", dict(buckets))
print("submissions by side", out["submissions_by_side"])
for k in ("other", "area", "vague"):
    print(f"\n== top {k} reported phenotypes")
    for t, n in side_terms[k].most_common(40):
        print(f"{n:6d}  {t[:150]}")
