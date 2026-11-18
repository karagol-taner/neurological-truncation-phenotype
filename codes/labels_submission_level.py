import csv, gzip, json, os, re, collections
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.dirname(HERE)
VS = os.path.join(SCRATCH, "variant_summary.txt.gz")
SS = os.path.join(SCRATCH, "submission_summary.txt.gz")
PH = json.load(open(os.path.join(HERE, "derived_phenotypes.json")))
GN = json.load(open(os.path.join(HERE, "gnomad_plof.json")))
V2 = json.load(open(os.path.join(HERE, "derived_v2.json")))

from area import IN_AREA
VAGUE = re.compile(r'not provided|not specified|see cases|none provided|not applicable|unknown|'
                   r'^-$|related disorder', re.I)
PTC = re.compile(r'p\.[A-Za-z]{3}(\d+)(Ter|\*)', re.I)
FS = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs', re.I)

POOL = sorted(g for g, r in PH["dual_phenotype_genes"].items() if r["positioned_ptv"] >= 10)
sym2gene = {s: g for g, e in V2["genes"].items() for s in e["clinvar_symbols"]}

def ptv_pos(name):
    m = PTC.search(name) or FS.search(name)
    return int(m.group(1)) if m else None


vid = {}
with gzip.open(VS, "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("Assembly") != "GRCh38": continue
        hits = {sym2gene[s] for s in (row.get("GeneSymbol") or "").split(";")
                if s in sym2gene and sym2gene[s] in POOL}
        if not hits: continue
        sig = (row.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig: continue
        pos = ptv_pos(row.get("Name") or "")
        if not pos: continue
        vid[row["VariationID"]] = (sorted(hits)[0], pos)
print("pathogenic truncating variants in pool genes:", len(vid))

counts = collections.defaultdict(collections.Counter)
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
        ph = r.get("ReportedPhenotypeInfo") or ""
        ph = re.sub(r'^[A-Z0-9:]+[:;]', '', ph)
        if VAGUE.search(ph) and not IN_AREA.search(ph):
            counts[v]["vague"] += 1
            continue
        counts[v]["area" if IN_AREA.search(ph) else "other"] += 1

side = collections.defaultdict(lambda: collections.defaultdict(list))
bucket_tot = collections.Counter()
for v, (g, pos) in vid.items():
    c = counts[v]
    if c["area"] and not c["other"]: b = "area"
    elif c["other"] and not c["area"]: b = "other"
    elif c["area"] and c["other"]: b = "both"
    else: b = "unlabelled"
    bucket_tot[b] += 1
    if b in ("area", "other"): side[g][b].append(pos)
print("per-variant labels from submissions:", dict(bucket_tot),
      "| resolved fraction:",
      round((bucket_tot['area'] + bucket_tot['other']) / max(sum(bucket_tot.values()), 1), 3))

out = {"label_buckets": dict(bucket_tot), "genes": {}}
pooled = collections.defaultdict(lambda: [0, 0, 0, 0])
DEFS = ("nmd_50nt_rule", "last_exon", "last_10pct", "last_25pct")
for g in POOL:
    gn = GN["genes"].get(g) or {}
    cds = gn.get("cds") or {}
    L = cds.get("protein_aa") or V2["genes"][g].get("protein_length")
    thrs = {"nmd_50nt_rule": cds.get("nmd_escape_from_codon"),
            "last_exon": cds.get("last_exon_start_codon"),
            "last_10pct": int(0.9 * L) + 1 if L else None,
            "last_25pct": int(0.75 * L) + 1 if L else None}
    dent, nond = side[g]["area"], side[g]["other"]
    tol = [v["protein_pos"] for v in gn.get("plof", [])
           if v.get("protein_pos") and v.get("lof") == "HC" and v.get("pass")]
    rec = {"literature_symbol": gn.get("literature_symbol"), "protein_aa": L,
           "n_area": len(dent), "n_other": len(nond),
           "n_gnomad_hc_pass": len(tol), "windows": {}}
    for name in DEFS:
        thr = thrs.get(name)
        if not thr: continue
        a = sum(1 for x in dent if x >= thr); b = len(dent) - a
        c = sum(1 for x in nond if x >= thr); d = len(nond) - c
        e = sum(1 for x in tol if x >= thr); f = len(tol) - e
        w = {"threshold_rel": round(thr / L, 3) if L else None,
             "area_in": a, "area_out": b, "other_in": c, "other_out": d,
             "gnomad_in": e, "gnomad_out": f}
        if dent and nond: w["p_area_vs_other"] = float(fisher_exact([[a, b], [c, d]])[1])
        if dent and tol: w["p_area_vs_gnomad"] = float(fisher_exact([[a, b], [e, f]])[1])
        rec["windows"][name] = w
        if dent and nond:
            p = pooled[name]; pooled[name] = [p[0]+a, p[1]+b, p[2]+c, p[3]+d]
    out["genes"][g] = rec
out["pooled"] = {n: {"area_in": a, "area_out": b, "other_in": c, "other_out": d,
                     "p": float(fisher_exact([[a, b], [c, d]])[1])}
                 for n, (a, b, c, d) in pooled.items()}
json.dump(out, open(os.path.join(HERE, "derived_labels_submission.json"), "w"),
          indent=1, sort_keys=True)

for name in DEFS:
    if name not in out["pooled"]: continue
    print(f"\n=== {name}")
    for g in sorted(POOL, key=lambda x: -(out["genes"][x]["n_area"] + out["genes"][x]["n_other"])):
        r = out["genes"][g]; w = r["windows"].get(name)
        if not w: continue
        p1, p2 = w.get("p_area_vs_other"), w.get("p_area_vs_gnomad")
        print(f"{(r.get('literature_symbol') or g):9s} thr={str(w['threshold_rel']):>6s} "
              f"D={w['area_in']}/{w['area_out']:<4d} N={w['other_in']}/{w['other_out']:<4d} "
              f"p(D/N)={(f'{p1:.2e}' if p1 is not None else 'na'):>9s} "
              f"G={w['gnomad_in']}/{w['gnomad_out']:<4d} p(D/G)={(f'{p2:.2e}' if p2 is not None else 'na'):>9s}")
    pl = out["pooled"][name]
    print(f"  pooled p = {pl['p']:.3e} (area {pl['area_in']}/{pl['area_out']}, "
          f"other {pl['other_in']}/{pl['other_out']})")
sig = [(out["genes"][g].get("literature_symbol") or g) for g in POOL
       if any((w.get("p_area_vs_other") or 1) < 0.05 for w in out["genes"][g]["windows"].values())]
print("\ngenes separating at p < 0.05 in at least one window definition:", len(sig), sig)
