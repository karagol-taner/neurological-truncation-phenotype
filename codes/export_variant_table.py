import csv, gzip, json, os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.dirname(HERE)
VS = os.path.join(SCRATCH, "variant_summary.txt.gz")
SS = os.path.join(SCRATCH, "submission_summary.txt.gz")
GN = json.load(open(os.path.join(HERE, "gnomad_plof.json")))
V2 = json.load(open(os.path.join(HERE, "derived_v2.json")))
TX = json.load(open(os.path.join(HERE, "tx_cache.json")))
OUT = os.path.join(HERE, "variant_table.csv")

GENES = set(g for g, e in V2["genes"].items() if e["tier"] == "isolated")
sym2gene = {s: g for g, e in V2["genes"].items() for s in e["clinvar_symbols"]}
LIT = {g: (GN["genes"].get(g) or {}).get("literature_symbol") or g for g in GENES}

from area import IN_AREA
VAGUE = re.compile(r'not provided|not specified|see cases|none provided|not applicable|unknown|'
                   r'^-$|related disorder', re.I)
PTC = re.compile(r'p\.[A-Za-z]{3}(\d+)(Ter|\*)', re.I)
FS = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs', re.I)
NAME = re.compile(r'^(N[MR]_\d+\.\d+)\(([^)]+)\):c\.([^ ]+)')
FSTER = re.compile(r'fsTer(\d+)', re.I)
DEL = re.compile(r'^(\d+)(?:_(\d+))?del([ACGT]*)$')
DUP = re.compile(r'^(\d+)(?:_(\d+))?dup([ACGT]*)$')
INS = re.compile(r'^(\d+)_(\d+)ins([ACGT]+)$')
DELINS = re.compile(r'^(\d+)(?:_(\d+))?delins([ACGT]+)$')
B = "ACGT"
CODON = {}
for i, aa in enumerate("KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"):
    CODON[B[i // 16] + B[(i // 4) % 4] + B[i % 4]] = aa
STOPS = {c for c, a in CODON.items() if a == "*"}

def apply_change(seq, c, cds_start):
    g = lambda p: cds_start + int(p) - 1
    m = DEL.match(c)
    if m:
        a, b = g(m.group(1)), g(m.group(2) or m.group(1)); return seq[:a-1] + seq[b:], a
    m = DUP.match(c)
    if m:
        a, b = g(m.group(1)), g(m.group(2) or m.group(1)); return seq[:b] + seq[a-1:b] + seq[b:], a
    m = INS.match(c)
    if m:
        a = g(m.group(1)); return seq[:a] + m.group(3) + seq[a:], a + 1
    m = DELINS.match(c)
    if m:
        a, b = g(m.group(1)), g(m.group(2) or m.group(1)); return seq[:a-1] + m.group(3) + seq[b:], a
    return None, None

def fsd(acc, c):
    rec = TX.get(acc)
    if not rec:
        return None
    mut, off = apply_change(rec["mrna"], c, rec["cds_start"])
    if mut is None or (len(mut) - len(rec["mrna"])) % 3 == 0:
        return None
    start = rec["cds_start"]
    i = start - 1 + ((off - start) // 3) * 3
    n = 0
    while i + 3 <= len(mut):
        cd = mut[i:i+3]
        if cd in STOPS:
            return n
        if cd not in CODON:
            return None
        n += 1; i += 3
    return None


var = {}
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
        m = FS.search(name)
        cls, pos = ("frameshift", int(m.group(1))) if m else (None, None)
        if not cls:
            m = PTC.search(name)
            cls, pos = ("nonsense", int(m.group(1))) if m else (None, None)
        if not cls or not pos:
            continue
        acc = chg = None
        mm = NAME.match(name)
        if mm:
            acc, chg = mm.group(1), mm.group(3).split()[0]
        var[row["VariationID"]] = {"gene": sorted(hits)[0], "cls": cls, "codon": pos,
                                   "acc": acc, "c": chg,
                                   "fsd": (fsd(acc, chg) if (cls == "frameshift" and acc) else None)}


cnt = collections.defaultdict(collections.Counter)
with gzip.open(SS, "rt", encoding="utf-8", errors="replace") as fh:
    header = None
    for line in fh:
        if line.startswith("#"):
            if "VariationID" in line:
                header = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if header is None:
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < len(header):
            continue
        r = dict(zip(header, f))
        v = r.get("VariationID")
        if v not in var:
            continue
        sig = (r.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig:
            continue
        ph = re.sub(r'^[A-Z0-9:]+[:;]', '', r.get("ReportedPhenotypeInfo") or "")
        if VAGUE.search(ph) and not IN_AREA.search(ph):
            continue
        cnt[v]["area" if IN_AREA.search(ph) else "other"] += 1

COLS = ["gene", "literature_symbol", "source", "variant_class", "codon", "protein_aa",
        "rel_pos", "phenotype_side", "fsd", "allele_count", "homozygote_count",
        "nmd_escape_from_codon", "nmd_escape_rel", "in_escape_window"]
rows_out = []
for vid, r in var.items():
    g = r["gene"]
    cds = (GN["genes"].get(g) or {}).get("cds") or {}
    L = cds.get("protein_aa") or V2["genes"][g].get("protein_length")
    esc = cds.get("nmd_escape_from_codon")
    c = cnt[vid]
    side = ("area" if c["area"] and not c["other"] else
            "other" if c["other"] and not c["area"] else
            "both" if c["area"] and c["other"] else "unlabelled")
    rows_out.append([g, LIT[g], "clinvar_pathogenic", r["cls"], r["codon"], L,
                     round(r["codon"] / L, 4) if L else "", side,
                     r["fsd"] if r["fsd"] is not None else "", "", "",
                     esc or "", round(esc / L, 4) if (esc and L) else "",
                     int(r["codon"] >= esc) if esc else ""])
for g, rec in GN["genes"].items():
    key = rec.get("clinvar_symbol") or g
    if key not in GENES:
        continue
    cds = rec.get("cds") or {}
    L = cds.get("protein_aa") or V2["genes"][key].get("protein_length")
    esc = cds.get("nmd_escape_from_codon")
    for v in rec.get("plof", []):
        if v.get("lof") != "HC" or not v.get("pass"):
            continue
        if v["consequence"] not in ("stop_gained", "frameshift_variant"):
            continue
        cls = "nonsense" if v["consequence"] == "stop_gained" else "frameshift"
        pos = v.get("protein_pos")
        m = FSTER.search(v.get("hgvsp") or "")
        rows_out.append([key, LIT[key], "gnomad_plof", cls, pos or "", L,
                         round(pos / L, 4) if (pos and L) else "", "population",
                         int(m.group(1)) - 1 if m else "", v["ac"], v["homozygote_count"],
                         esc or "", round(esc / L, 4) if (esc and L) else "",
                         int(pos >= esc) if (pos and esc) else ""])
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(COLS); w.writerows(rows_out)

n = collections.Counter(r[2] for r in rows_out)
sides = collections.Counter(r[7] for r in rows_out if r[2] == "clinvar_pathogenic")
print("wrote", OUT, "rows:", len(rows_out))
print("by source:", dict(n))
print("pathogenic by phenotype side:", dict(sides))
print("with an FSD value:", sum(1 for r in rows_out if r[8] != ""))
