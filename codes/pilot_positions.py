import csv, gzip, json, os, re, collections
from scipy.stats import mannwhitneyu

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
    if m:
        return "nonsense", int(m.group(1))
    m = FS.search(name)
    if m:
        return "frameshift", int(m.group(1))
    return None, None

side = collections.defaultdict(lambda: collections.defaultdict(list))
with gzip.open(VS, "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("Assembly") != "GRCh38":
            continue
        hits = {sym2gene[s] for s in (row.get("GeneSymbol") or "").split(";")
                if s in sym2gene and sym2gene[s] in POOL}
        if not hits:
            continue
        sig = (row.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig:
            continue
        cls, pos = ptv_pos(row.get("Name") or "")
        if not pos:
            continue
        terms = [t.strip() for t in (row.get("PhenotypeList") or "").split("|") if t.strip()]
        useful = [t for t in terms if not VAGUE.match(t)]
        d = [t for t in useful if IN_AREA.search(t)]
        nd = [t for t in useful if not IN_AREA.search(t)]
        if d and not nd:
            b = "area"
        elif nd and not d:
            b = "other"
        else:
            continue
        for g in hits:
            side[g][b].append((pos, cls))

def med(xs):
    xs = sorted(xs)
    return None if not xs else round(xs[len(xs) // 2], 3)

out = {"pool": POOL, "genes": {}}
for g in POOL:
    gn = GN["genes"].get(g) or {}
    cds = gn.get("cds") or {}
    L = cds.get("protein_aa") or V2["genes"][g].get("protein_length")
    esc = cds.get("nmd_escape_from_codon")
    rec = {"protein_aa": L, "length_source": "MANE CDS" if cds.get("protein_aa") else "UniProt",
           "nmd_escape_from_codon": esc,
           "nmd_escape_rel": round(esc / L, 3) if (esc and L) else None,
           "literature_symbol": gn.get("literature_symbol")}
    sets = {}
    for b in ("area", "other"):
        p = [x for x, _ in side[g][b]]
        sets[b] = p
        rec[f"n_{b}"] = len(p)
        rec[f"rel_median_{b}"] = med([x / L for x in p]) if (p and L) else None
        rec[f"escape_frac_{b}"] = (round(sum(1 for x in p if x >= esc) / len(p), 3)
                                   if (p and esc) else None)
    tol = [v["protein_pos"] for v in gn.get("plof", [])
           if v.get("protein_pos") and v.get("lof") == "HC" and v.get("pass")]
    tol2 = [v["protein_pos"] for v in gn.get("plof", [])
            if v.get("protein_pos") and v.get("lof") == "HC" and v.get("pass") and v["ac"] >= 2]
    hom = sum(1 for v in gn.get("plof", []) if v.get("homozygote_count", 0) > 0)
    rec.update({"n_gnomad_hc_pass": len(tol), "n_gnomad_ac_ge2": len(tol2),
                "n_gnomad_with_homozygote": hom,
                "rel_median_gnomad": med([x / L for x in tol]) if (tol and L) else None,
                "escape_frac_gnomad": (round(sum(1 for x in tol if x >= esc) / len(tol), 3)
                                       if (tol and esc) else None)})
    if len(sets["area"]) >= 3 and len(sets["other"]) >= 3:
        u, p = mannwhitneyu(sets["area"], sets["other"], alternative="two-sided")
        rec["mwu_p_area_vs_other"] = float(p)
    if len(sets["area"]) >= 3 and len(tol) >= 3:
        rec["mwu_p_area_vs_gnomad"] = float(
            mannwhitneyu(sets["area"], tol, alternative="two-sided")[1])
    out["genes"][g] = rec

json.dump(out, open(os.path.join(HERE, "derived_pilot.json"), "w"), indent=1, sort_keys=True)

print(f"{'gene':9s} {'aa':>5s} {'escRel':>7s} {'nDent':>6s} {'relD':>6s} {'escD':>5s} "
      f"{'nNonD':>6s} {'relN':>6s} {'escN':>5s} {'p(D/N)':>9s} {'gnomHC':>7s} {'relG':>6s} {'escG':>5s}")
for g in sorted(POOL, key=lambda x: -(out["genes"][x]["n_area"] + out["genes"][x]["n_other"])):
    r = out["genes"][g]
    name = r.get("literature_symbol") or g
    pv = r.get("mwu_p_area_vs_other")
    print(f"{name:9s} {str(r['protein_aa']):>5s} {str(r['nmd_escape_rel']):>7s} "
          f"{r['n_area']:6d} {str(r['rel_median_area']):>6s} {str(r['escape_frac_area']):>5s} "
          f"{r['n_other']:6d} {str(r['rel_median_other']):>6s} {str(r['escape_frac_other']):>5s} "
          f"{(f'{pv:.2e}' if pv is not None else 'na'):>9s} "
          f"{r['n_gnomad_hc_pass']:7d} {str(r['rel_median_gnomad']):>6s} {str(r['escape_frac_gnomad']):>5s}")
sig = [g for g in POOL if (out["genes"][g].get("mwu_p_area_vs_other") or 1) < 0.05]
print(f"\ngenes with a two-sided separation at p < 0.05: {len(sig)} of {len(POOL)}: "
      f"{', '.join(out['genes'][g].get('literature_symbol') or g for g in sig)}")
