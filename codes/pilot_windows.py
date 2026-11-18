import csv, gzip, json, os, re, collections
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.dirname(HERE)
VS = os.path.join(SCRATCH, "variant_summary.txt.gz")
PH = json.load(open(os.path.join(HERE, "derived_phenotypes.json")))
GN = json.load(open(os.path.join(HERE, "gnomad_plof.json")))
V2 = json.load(open(os.path.join(HERE, "derived_v2.json")))

from area import IN_AREA
VAGUE = re.compile(r'^(not provided|not specified|see cases|none provided|not applicable|unknown)$', re.I)
PTC = re.compile(r'p\.[A-Za-z]{3}(\d+)(Ter|\*)', re.I)
FS = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs', re.I)

POOL = sorted(g for g, r in PH["dual_phenotype_genes"].items() if r["positioned_ptv"] >= 10)
sym2gene = {s: g for g, e in V2["genes"].items() for s in e["clinvar_symbols"]}

def ptv_pos(name):
    m = PTC.search(name)
    if m: return int(m.group(1))
    m = FS.search(name)
    if m: return int(m.group(1))
    return None

side = collections.defaultdict(lambda: collections.defaultdict(list))
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
        useful = [t.strip() for t in (row.get("PhenotypeList") or "").split("|")
                  if t.strip() and not VAGUE.match(t.strip())]
        d = [t for t in useful if IN_AREA.search(t)]
        nd = [t for t in useful if not IN_AREA.search(t)]
        b = "area" if (d and not nd) else ("other" if (nd and not d) else None)
        if not b: continue
        for g in hits: side[g][b].append(pos)

def windows(cds, L):

    w = {}
    if cds.get("nmd_escape_from_codon"): w["nmd_50nt_rule"] = cds["nmd_escape_from_codon"]
    if cds.get("last_exon_start_codon"): w["last_exon"] = cds["last_exon_start_codon"]
    if L:
        w["last_10pct"] = int(0.90 * L) + 1
        w["last_25pct"] = int(0.75 * L) + 1
    return w

out = {"definitions": ["nmd_50nt_rule", "last_exon", "last_10pct", "last_25pct"], "genes": {}}
pooled = collections.defaultdict(lambda: [0, 0, 0, 0])
for g in POOL:
    gn = GN["genes"].get(g) or {}
    cds = gn.get("cds") or {}
    L = cds.get("protein_aa") or V2["genes"][g].get("protein_length")
    dent, nond = side[g]["area"], side[g]["other"]
    tol = [v["protein_pos"] for v in gn.get("plof", [])
           if v.get("protein_pos") and v.get("lof") == "HC" and v.get("pass")]
    rec = {"literature_symbol": gn.get("literature_symbol"), "protein_aa": L,
           "n_area": len(dent), "n_other": len(nond), "n_gnomad_hc_pass": len(tol),
           "windows": {}}
    for name, thr in windows(cds, L).items():
        a = sum(1 for x in dent if x >= thr); b = len(dent) - a
        c = sum(1 for x in nond if x >= thr); d = len(nond) - c
        e = sum(1 for x in tol if x >= thr); f = len(tol) - e
        w = {"threshold_codon": thr, "threshold_rel": round(thr / L, 3) if L else None,
             "area_in": a, "area_out": b, "other_in": c, "other_out": d,
             "gnomad_in": e, "gnomad_out": f}
        if dent and nond:
            w["p_area_vs_other"] = float(fisher_exact([[a, b], [c, d]])[1])
        if dent and tol:
            w["p_area_vs_gnomad"] = float(fisher_exact([[a, b], [e, f]])[1])
        rec["windows"][name] = w
        if dent and nond:
            p = pooled[name]
            pooled[name] = [p[0] + a, p[1] + b, p[2] + c, p[3] + d]
    out["genes"][g] = rec

out["pooled"] = {}
for name, (a, b, c, d) in pooled.items():
    out["pooled"][name] = {"area_in": a, "area_out": b, "other_in": c,
                           "other_out": d,
                           "p": float(fisher_exact([[a, b], [c, d]])[1])}
json.dump(out, open(os.path.join(HERE, "derived_pilot_windows.json"), "w"), indent=1, sort_keys=True)

for name in out["definitions"]:
    print(f"\n=== window definition: {name}")
    print(f"{'gene':9s} {'thrRel':>7s} {'dent in/out':>12s} {'nonD in/out':>12s} "
          f"{'p(D/N)':>10s} {'gnomAD in/out':>14s} {'p(D/gnomAD)':>12s}")
    for g in sorted(POOL, key=lambda x: -(out["genes"][x]["n_area"] + out["genes"][x]["n_other"])):
        r = out["genes"][g]; w = r["windows"].get(name)
        if not w: continue
        nm = r.get("literature_symbol") or g
        p1 = w.get("p_area_vs_other"); p2 = w.get("p_area_vs_gnomad")
        print(f"{nm:9s} {str(w['threshold_rel']):>7s} "
              f"{str(w['area_in'])+'/'+str(w['area_out']):>12s} "
              f"{str(w['other_in'])+'/'+str(w['other_out']):>12s} "
              f"{(f'{p1:.2e}' if p1 is not None else 'na'):>10s} "
              f"{str(w['gnomad_in'])+'/'+str(w['gnomad_out']):>14s} "
              f"{(f'{p2:.2e}' if p2 is not None else 'na'):>12s}")
    pl = out["pooled"][name]
    print(f"  pooled: area {pl['area_in']}/{pl['area_out']}, "
          f"other {pl['other_in']}/{pl['other_out']}, p = {pl['p']:.3e}")
