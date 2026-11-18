import csv, gzip, json, os, re, statistics, sys, time, urllib.request, collections
from scipy.stats import mannwhitneyu, fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.dirname(HERE)
VS = os.path.join(SCRATCH, "variant_summary.txt.gz")
CACHE = os.path.join(HERE, "tx_cache.json")
GN = json.load(open(os.path.join(HERE, "gnomad_plof.json")))
V2 = json.load(open(os.path.join(HERE, "derived_v2.json")))
GENES = set(g for g, e in V2["genes"].items() if e["tier"] == "isolated")
sym2gene = {s: g for g, e in V2["genes"].items() for s in e["clinvar_symbols"]}

NAME = re.compile(r'^(N[MR]_\d+\.\d+)\(([^)]+)\):c\.([^ ]+)')
DEL = re.compile(r'^(\d+)(?:_(\d+))?del([ACGT]*)$')
DUP = re.compile(r'^(\d+)(?:_(\d+))?dup([ACGT]*)$')
INS = re.compile(r'^(\d+)_(\d+)ins([ACGT]+)$')
DELINS = re.compile(r'^(\d+)(?:_(\d+))?delins([ACGT]+)$')
FSTER = re.compile(r'fsTer(\d+)', re.I)
B = "ACGT"
CODON = {}
for i, aa in enumerate("KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"):
    CODON[B[i // 16] + B[(i // 4) % 4] + B[i % 4]] = aa
STOPS = {c for c, a in CODON.items() if a == "*"}

cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}

def fetch(acc):

    if acc in cache:
        return cache[acc]
    def get(rettype):
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="
               + acc + "&rettype=" + rettype + "&retmode=text")
        for a in range(4):
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    return r.read().decode()
            except Exception as e:
                if a == 3:
                    print("fetch failed", acc, rettype, e, file=sys.stderr)
                    return None
                time.sleep(3 * (a + 1))
    fa, cdsfa = get("fasta"), get("fasta_cds_na")
    rec = None
    if fa and cdsfa:
        m = "".join(l.strip() for l in fa.splitlines() if not l.startswith(">")).upper()
        c = "".join(l.strip() for l in cdsfa.splitlines() if not l.startswith(">")).upper()
        i = m.find(c)
        if i >= 0 and c:
            rec = {"mrna": m, "cds_start": i + 1, "cds_len": len(c)}
    cache[acc] = rec
    json.dump(cache, open(CACHE, "w"))
    time.sleep(0.12)
    return rec

def apply_change(seq, c, cds_start):

    def g(p):
        return cds_start + int(p) - 1
    m = DEL.match(c)
    if m:
        a, b = g(m.group(1)), g(m.group(2) or m.group(1))
        return seq[:a - 1] + seq[b:], a
    m = DUP.match(c)
    if m:
        a, b = g(m.group(1)), g(m.group(2) or m.group(1))
        return seq[:b] + seq[a - 1:b] + seq[b:], a
    m = INS.match(c)
    if m:
        a = g(m.group(1))
        return seq[:a] + m.group(3) + seq[a:], a + 1
    m = DELINS.match(c)
    if m:
        a, b = g(m.group(1)), g(m.group(2) or m.group(1))
        return seq[:a - 1] + m.group(3) + seq[b:], a
    return None, None

def fsd(rec, c):

    mut, off = apply_change(rec["mrna"], c, rec["cds_start"])
    if mut is None:
        return None, None
    if (len(mut) - len(rec["mrna"])) % 3 == 0:
        return None, None
    start = rec["cds_start"]
    codon_idx = (off - start) // 3
    i = start - 1 + codon_idx * 3
    novel = 0
    while i + 3 <= len(mut):
        cd = mut[i:i + 3]
        if cd in STOPS:
            return novel, False
        if cd not in CODON:
            return None, None
        novel += 1
        i += 3
    return novel, True

path = collections.defaultdict(list)
skipped = collections.Counter()
seen = set()
with gzip.open(VS, "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("Assembly") != "GRCh38":
            continue
        hits = {sym2gene[s] for s in (row.get("GeneSymbol") or "").split(";")
                if s in sym2gene and sym2gene[s] in GENES}
        if not hits:
            continue
        sig = (row.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig:
            continue
        name = row.get("Name") or ""
        if "fs" not in name:
            continue
        m = NAME.match(name)
        if not m:
            skipped["name_unparsed"] += 1
            continue
        acc, chg = m.group(1), m.group(3).split()[0]
        if (acc, chg) in seen:
            continue
        seen.add((acc, chg))
        rec = fetch(acc)
        if not rec:
            skipped["no_transcript"] += 1
            continue
        novel, cens = fsd(rec, chg)
        if novel is None:
            skipped["change_unparsed_or_inframe"] += 1
            continue
        for g in hits:
            path[g].append({"fsd": novel, "censored": bool(cens), "c": chg, "acc": acc})

pop = collections.defaultdict(list)
for g, rec in GN["genes"].items():
    key = rec.get("clinvar_symbol") or g
    if key not in GENES:
        continue
    for v in rec.get("plof", []):
        if v["consequence"] != "frameshift_variant" or v.get("lof") != "HC" or not v.get("pass"):
            continue
        m = FSTER.search(v.get("hgvsp") or "")
        if m:
            pop[key].append({"fsd": int(m.group(1)) - 1, "censored": False, "ac": v["ac"]})

THRESH = (10, 20, 40, 80)
out = {"definition": "novel amino acids translated before the new termination codon",
       "note": "pathogenic side computed on the full transcript so a neo-tail may run past the "
               "original stop; population side read from gnomAD's own protein notation",
       "skipped": dict(skipped), "thresholds": list(THRESH), "genes": {}}
allp, allq, cens_tot = [], [], 0
for g in sorted(GENES):
    P = [r["fsd"] for r in path[g]]
    C = sum(1 for r in path[g] if r["censored"])
    Q = [r["fsd"] for r in pop[g]]
    cens_tot += C
    med = lambda xs: (statistics.median(xs) if xs else None)
    rec = {"literature_symbol": (GN["genes"].get(g) or {}).get("literature_symbol"),
           "n_pathogenic": len(P), "n_pathogenic_censored": C, "n_population": len(Q),
           "median_fsd_pathogenic": med(P), "median_fsd_population": med(Q),
           "max_fsd_pathogenic": max(P) if P else None,
           "max_fsd_population": max(Q) if Q else None, "long": {}}
    if len(P) >= 3 and len(Q) >= 3:
        rec["mwu_p"] = float(mannwhitneyu(P, Q, alternative="two-sided")[1])
        for t in THRESH:
            a = sum(1 for x in P if x >= t); b = len(P) - a
            c = sum(1 for x in Q if x >= t); d = len(Q) - c
            rec["long"][f"ge{t}"] = {"pathogenic": f"{a}/{len(P)}", "population": f"{c}/{len(Q)}",
                                     "p": float(fisher_exact([[a, b], [c, d]])[1])}
    out["genes"][g] = rec
    allp += P; allq += Q
out["pooled"] = {"n_pathogenic": len(allp), "n_pathogenic_censored": cens_tot,
                 "n_population": len(allq),
                 "median_pathogenic": statistics.median(allp) if allp else None,
                 "median_population": statistics.median(allq) if allq else None}
if len(allp) >= 3 and len(allq) >= 3:
    out["pooled"]["mwu_p"] = float(mannwhitneyu(allp, allq, alternative="two-sided")[1])
    for t in THRESH:
        a = sum(1 for x in allp if x >= t); b = len(allp) - a
        c = sum(1 for x in allq if x >= t); d = len(allq) - c
        out["pooled"][f"long_ge{t}"] = {"pathogenic": f"{a}/{len(allp)}",
                                        "population": f"{c}/{len(allq)}",
                                        "p": float(fisher_exact([[a, b], [c, d]])[1])}
json.dump(out, open(os.path.join(HERE, "derived_fsd2.json"), "w"), indent=1, sort_keys=True)

print("skipped:", dict(skipped), "| right-censored pathogenic:", cens_tot)
print(f"{'gene':9s} {'nP':>4s} {'cens':>5s} {'nQ':>4s} {'medP':>5s} {'medQ':>5s} {'maxP':>5s} "
      f"{'MWU p':>9s} " + " ".join(f"{'p>='+str(t):>9s}" for t in THRESH))
for g in sorted(out["genes"], key=lambda x: -out["genes"][x]["n_pathogenic"]):
    r = out["genes"][g]
    if not r["n_pathogenic"]:
        continue
    nm = r.get("literature_symbol") or g
    pv = f"{r['mwu_p']:.2e}" if r.get("mwu_p") is not None else "na"
    longs = " ".join((f"{r['long'][f'ge{t}']['p']:.2e}" if f"ge{t}" in r["long"] else "na").rjust(9)
                     for t in THRESH)
    print(f"{nm:9s} {r['n_pathogenic']:4d} {r['n_pathogenic_censored']:5d} {r['n_population']:4d} "
          f"{str(r['median_fsd_pathogenic']):>5s} {str(r['median_fsd_population']):>5s} "
          f"{str(r['max_fsd_pathogenic']):>5s} {pv:>9s} {longs}")
print("\npooled:", json.dumps(out["pooled"]))
