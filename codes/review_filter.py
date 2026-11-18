import csv, gzip, json, os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
VS = os.path.join(os.path.dirname(HERE), "variant_summary.txt.gz")
V2 = json.load(open(os.path.join(HERE, "derived_v2.json")))

canon = {}
for g, e in V2["genes"].items():
    for s in e["clinvar_symbols"]:
        canon[s] = g
want = set(canon)

PTC = re.compile(r'p\.[A-Za-z]{3}(\d+)(Ter|\*)', re.I)
FS = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs', re.I)

def ptv(name):
    m = PTC.search(name)
    if m:
        return "nonsense", int(m.group(1))
    m = FS.search(name)
    if m:
        return "frameshift", int(m.group(1))
    return None, None

LEVELS = ("all", "criteria", "multi_or_expert")

def levels_of(review, nsub):
    out = ["all"]
    r = (review or "").lower()
    if "criteria provided" in r or "expert panel" in r or "practice guideline" in r:
        out.append("criteria")
    if "expert panel" in r or "practice guideline" in r or "multiple submitters" in r or nsub >= 2:
        out.append("multi_or_expert")
    return out

counts = collections.defaultdict(collections.Counter)
pos = collections.defaultdict(lambda: collections.defaultdict(list))

with gzip.open(VS, "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("Assembly") != "GRCh38":
            continue
        hits = {canon[s] for s in (row.get("GeneSymbol") or "").split(";") if s in want}
        if not hits:
            continue
        sig = (row.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig:
            continue
        cls, p = ptv(row.get("Name") or "")
        if not cls:
            continue
        try:
            nsub = int(row.get("NumberSubmitters") or 0)
        except ValueError:
            nsub = 0
        lv = levels_of(row.get("ReviewStatus"), nsub)
        for g in hits:
            for l in lv:
                counts[g][l] += 1
                counts[g][f"{l}:{cls}"] += 1
                if p:
                    pos[g][l].append(p)

out = {"levels": list(LEVELS), "genes": {}}
for g, e in V2["genes"].items():
    L = e.get("protein_length")
    rec = {"tier": e["tier"], "protein_length": L}
    for l in LEVELS:
        pp = sorted(pos[g][l])
        rec[f"{l}_ptv"] = counts[g][l]
        rec[f"{l}_nonsense"] = counts[g][f"{l}:nonsense"]
        rec[f"{l}_frameshift"] = counts[g][f"{l}:frameshift"]
        rec[f"{l}_positioned"] = len(pp)
        if pp and L:
            rel = sorted(x / L for x in pp)
            rec[f"{l}_rel_median"] = round(rel[len(rel) // 2], 3)
            rec[f"{l}_frac_gt90"] = round(sum(1 for x in rel if x > 0.9) / len(rel), 3)
    out["genes"][g] = rec

iso = [r for r in out["genes"].values() if r["tier"] == "isolated"]
for l in LEVELS:
    out[f"isolated_total_{l}"] = sum(r[f"{l}_ptv"] for r in iso)
    for n in (5, 10, 20):
        out[f"isolated_genes_{l}_ge{n}"] = sum(1 for r in iso if r[f"{l}_positioned"] >= n)
json.dump(out, open(os.path.join(HERE, "derived_review.json"), "w"), indent=1, sort_keys=True)

print("isolated PTV totals by evidence level:")
for l in LEVELS:
    print(f"  {l:16s} total={out[f'isolated_total_{l}']:5d} "
          f">=5:{out[f'isolated_genes_{l}_ge5']:3d} "
          f">=10:{out[f'isolated_genes_{l}_ge10']:3d} "
          f">=20:{out[f'isolated_genes_{l}_ge20']:3d}")
