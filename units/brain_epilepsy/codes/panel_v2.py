import csv, gzip, json, os, re, sys, collections, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.dirname(HERE)
GCS = os.path.join(SCRATCH, "gene_condition_source_id")
VS = os.path.join(SCRATCH, "variant_summary.txt.gz")

DISEASE = re.compile('epilep|seizure|epileptic encephalopathy|infantile spasms|dravet|lennox|ohtahara|west syndrome', re.I)
EXCL = re.compile('neuroblastoma|glioma|glioblastoma|medulloblastoma|meningioma|ependymoma|astrocytoma|cerebrotendinous|cerebro-oculo-facio|cerebrocostomandibular|ciliary dyskinesia', re.I)
SYNDROMIC = re.compile('dysmorph|facial|facies|short stature|growth (retardation|delay)|polydactyly|syndactyly|skeletal|limb (anomal|defect|malform|reduction|deficiency)|obesity|immunodeficien|deafness|hearing|retin|optic atrophy|cataract|ocular|ophthalm|renal|kidney|nephr|cardi|heart|hepat|liver|ectodermal|hair|albinism|myopathy|muscular dystrophy|neuropathy|hypogonad|diabetes|anemia|thrombocytopenia|cleft|craniosynostosis|osteo|chondro|lung cysts|wilms tumor', re.I)
ISOLATED = re.compile('^(developmental and epileptic encephalopathy|epileptic encephalopathy|early infantile epileptic|epilepsy|generalized epilepsy|focal epilepsy|familial (focal|temporal|adult myoclonic|partial|infantile|benign)|benign familial (neonatal|infantile)|juvenile myoclonic|childhood absence|febrile seizures|seizures|myoclonic epilepsy|progressive myoclonic epilepsy|generalized epilepsy with febrile|infantile spasms)', re.I)

PTC = re.compile(r'p\.[A-Za-z]{3}(\d+)(Ter|\*)', re.I)
FS = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs', re.I)
MIS = re.compile(r'p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})')
SPLICE = re.compile(r'c\.[-*]?\d+[+-]\d+')


def panel():
    genes = {}
    with open(GCS, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            dis = (r.get("DiseaseName") or "").strip()
            sym = (r.get("AssociatedGenes") or "").strip()
            if not sym or not DISEASE.search(dis) or EXCL.search(dis):
                continue
            rec = genes.setdefault(sym, {"conditions": set()})
            rec["conditions"].add(dis)
    for sym, rec in genes.items():
        iso = [c for c in rec["conditions"]
               if ISOLATED.match(c) and not SYNDROMIC.search(c)]
        rec["isolated"] = sorted(iso)
        rec["syndromic"] = sorted(c for c in rec["conditions"] if c not in set(iso))
        rec["tier"] = "isolated" if iso else "syndromic_only"
    return genes


def uniprot(symbols):

    out = {}
    for s in symbols:
        q = f"gene:{s} AND organism_id:9606 AND reviewed:true"
        url = ("https://rest.uniprot.org/uniprotkb/search?query="
               + urllib.parse.quote(q)
               + "&fields=accession,length,gene_primary&format=tsv&size=5")
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                rows = r.read().decode().splitlines()[1:]
        except Exception as e:
            print("uniprot fail", s, e, file=sys.stderr)
            continue
        if not rows:
            continue
        f = rows[0].split("\t")
        if len(f) >= 3:
            try:
                out[s] = {"acc": f[0], "length": int(f[1]), "primary": f[2].strip()}
            except ValueError:
                pass
    return out


def classify(name, vtype):
    m = PTC.search(name)
    if m:
        return "nonsense", int(m.group(1))
    m = FS.search(name)
    if m:
        return "frameshift", int(m.group(1))
    if SPLICE.search(name) and "p." not in name:
        return "splice_region", None
    m = MIS.search(name)
    if m:
        a, pos, b = m.group(1), int(m.group(2)), m.group(3)
        if a.lower() == b.lower():
            return "synonymous", pos
        if b.lower() == "ter":
            return "nonsense", pos
        return "missense", pos
    return "other", None


def bucket(s):
    s = (s or "").lower()
    if "conflicting" in s:
        return "conflicting"
    if "uncertain" in s:
        return "vus"
    if "pathogenic" in s:
        return "plp"
    if "benign" in s:
        return "benign"
    return "other"


def main():
    pan = panel()
    up = uniprot(sorted(pan))

    by_acc = collections.defaultdict(list)
    for s, rec in up.items():
        by_acc[rec["acc"]].append(s)
    alias_merges = {}
    canon = {}
    for acc, syms in by_acc.items():
        prim = up[syms[0]]["primary"]
        keep = prim if prim in syms else sorted(syms)[0]
        for s in syms:
            canon[s] = keep
        if len(syms) > 1:
            alias_merges[keep] = sorted(s for s in syms if s != keep)
    for s in pan:
        canon.setdefault(s, s)

    counts = collections.defaultdict(collections.Counter)
    pos = collections.defaultdict(lambda: collections.defaultdict(list))
    n_rows = 0
    want = set(pan)
    with gzip.open(VS, "rt", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("Assembly") != "GRCh38":
                continue
            hits = {canon[s] for s in (row.get("GeneSymbol") or "").split(";")
                    if s in want}
            if not hits:
                continue
            n_rows += 1
            b = bucket(row.get("ClinicalSignificance"))
            cls, p = classify(row.get("Name") or "", row.get("Type") or "")
            for g in hits:
                counts[g][b] += 1
                counts[g][f"{b}:{cls}"] += 1
                if cls in ("nonsense", "frameshift") and p:
                    pos[g][b].append((p, cls))

    genes = {}
    for s, rec in pan.items():
        g = canon[s]
        e = genes.setdefault(g, {"gene": g, "clinvar_symbols": set(),
                                 "isolated": set(), "syndromic": set()})
        e["clinvar_symbols"].add(s)
        e["isolated"] |= set(rec["isolated"])
        e["syndromic"] |= set(rec["syndromic"])
    for g, e in genes.items():
        c, p = counts[g], pos[g]
        L = up.get(g, {}).get("length")
        e["clinvar_symbols"] = sorted(e["clinvar_symbols"])
        e["isolated"] = sorted(e["isolated"])
        e["syndromic"] = sorted(e["syndromic"])
        e["tier"] = "isolated" if e["isolated"] else "syndromic_only"
        e["uniprot"] = up.get(g, {}).get("acc")
        e["protein_length"] = L
        e["plp_total"] = c["plp"]
        e["plp_nonsense"] = c["plp:nonsense"]
        e["plp_frameshift"] = c["plp:frameshift"]
        e["plp_missense"] = c["plp:missense"]
        e["plp_splice_region"] = c["plp:splice_region"]
        e["plp_ptv"] = c["plp:nonsense"] + c["plp:frameshift"]
        e["vus_total"] = c["vus"]
        e["vus_missense"] = c["vus:missense"]
        e["benign_total"] = c["benign"]
        e["conflicting_total"] = c["conflicting"]
        pp = sorted(x for x, _ in p["plp"])
        e["plp_ptv_positioned"] = len(pp)
        e["vus_ptv_positioned"] = len(p["vus"])
        e["benign_ptv_positioned"] = len(p["benign"])
        if pp and L:
            rel = [x / L for x in pp]
            rel.sort()
            e["plp_ptv_rel_median"] = round(rel[len(rel) // 2], 3)
            e["plp_ptv_rel_min"] = round(rel[0], 3)
            e["plp_ptv_rel_max"] = round(rel[-1], 3)
            e["plp_ptv_frac_last_10pct"] = round(sum(1 for x in rel if x > 0.9) / len(rel), 3)
            e["plp_ptv_frac_last_25pct"] = round(sum(1 for x in rel if x > 0.75) / len(rel), 3)
            e["plp_ptv_span_rel"] = round(rel[-1] - rel[0], 3)
        if e["plp_ptv"]:
            e["plp_fs_fraction"] = round(e["plp_frameshift"] / e["plp_ptv"], 3)

    iso = [e for e in genes.values() if e["tier"] == "isolated"]
    out = {
        "clinvar_release_file": "variant_summary.txt.gz",
        "clinvar_rows_matched_grch38": n_rows,
        "panel_symbols_matched": len(pan),
        "panel_genes_after_alias_merge": len(genes),
        "alias_merges": alias_merges,
        "n_isolated": len(iso),
        "n_syndromic_only": len(genes) - len(iso),
        "isolated_plp_total": sum(e["plp_total"] for e in iso),
        "isolated_ptv_total": sum(e["plp_ptv"] for e in iso),
        "isolated_vus_total": sum(e["vus_total"] for e in iso),
        "isolated_genes_ptv_ge_20": sum(1 for e in iso if e["plp_ptv_positioned"] >= 20),
        "isolated_genes_ptv_ge_10": sum(1 for e in iso if e["plp_ptv_positioned"] >= 10),
        "isolated_genes_ptv_ge_5": sum(1 for e in iso if e["plp_ptv_positioned"] >= 5),
        "genes": genes,
    }
    json.dump(out, open(os.path.join(HERE, "derived_v2.json"), "w"),
              indent=1, sort_keys=True, default=list)

    cols = ["gene", "tier", "clinvar_symbols", "uniprot", "protein_length",
            "plp_total", "plp_ptv", "plp_nonsense", "plp_frameshift",
            "plp_fs_fraction", "plp_missense", "plp_ptv_positioned",
            "plp_ptv_rel_min", "plp_ptv_rel_median", "plp_ptv_rel_max",
            "plp_ptv_frac_last_10pct", "plp_ptv_frac_last_25pct",
            "vus_total", "vus_missense", "benign_total", "conflicting_total",
            "isolated", "syndromic"]
    with open(os.path.join(HERE, "panel_v2.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for g in sorted(genes, key=lambda x: -genes[x]["plp_ptv"]):
            e = genes[g]
            w.writerow(["; ".join(e[c]) if isinstance(e.get(c), list) else e.get(c, "")
                        for c in cols])

    print(f"symbols matched {len(pan)} -> genes {len(genes)} "
          f"(isolated {len(iso)}, syndromic-only {out['n_syndromic_only']}); "
          f"GRCh38 rows {n_rows}")
    print("alias merges:", alias_merges)
    print(f"isolated panel: P/LP {out['isolated_plp_total']}, PTV {out['isolated_ptv_total']}, "
          f"VUS {out['isolated_vus_total']}; genes with >=20 positioned PTV: "
          f"{out['isolated_genes_ptv_ge_20']}, >=10: {out['isolated_genes_ptv_ge_10']}, "
          f">=5: {out['isolated_genes_ptv_ge_5']}")
    print("\nISOLATED GENES, by PTV count")
    for e in sorted(iso, key=lambda r: -r["plp_ptv"]):
        print(f"{e['gene']:9s} L={str(e['protein_length']):>5s} P/LP={e['plp_total']:4d} "
              f"PTV={e['plp_ptv']:4d} fs%={e.get('plp_fs_fraction')} "
              f"relmed={e.get('plp_ptv_rel_median')} last10={e.get('plp_ptv_frac_last_10pct')} "
              f"last25={e.get('plp_ptv_frac_last_25pct')} VUS={e['vus_total']:5d}")


if __name__ == "__main__":
    main()
