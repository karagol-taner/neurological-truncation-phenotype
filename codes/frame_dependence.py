import csv, gzip, json, os, re, collections
from scipy.stats import fisher_exact, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.dirname(HERE)
VS = os.path.join(SCRATCH, "variant_summary.txt.gz")
GN = json.load(open(os.path.join(HERE, "gnomad_plof.json")))
V2 = json.load(open(os.path.join(HERE, "derived_v2.json")))

PTC = re.compile(r'p\.[A-Za-z]{3}(\d+)(Ter|\*)', re.I)
FS = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs', re.I)
GENES = sorted(g for g, e in V2["genes"].items() if e["tier"] == "isolated")
sym2gene = {s: g for g, e in V2["genes"].items() for s in e["clinvar_symbols"]}

def cls_pos(name):
    m = FS.search(name)
    if m:
        return "frameshift", int(m.group(1))
    m = PTC.search(name)
    if m:
        return "nonsense", int(m.group(1))
    return None, None

path = collections.defaultdict(lambda: {"nonsense": [], "frameshift": []})
with gzip.open(VS, "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("Assembly") != "GRCh38":
            continue
        hits = {sym2gene[s] for s in (row.get("GeneSymbol") or "").split(";")
                if s in sym2gene and sym2gene[s] in set(GENES)}
        if not hits:
            continue
        sig = (row.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig:
            continue
        c, p = cls_pos(row.get("Name") or "")
        if not c or not p:
            continue
        for g in hits:
            path[g][c].append(p)

def bh(pairs):

    xs = sorted(pairs, key=lambda kv: kv[1])
    n = len(xs)
    q, prev = {}, 1.0
    for i in range(n - 1, -1, -1):
        k, p = xs[i]
        prev = min(prev, p * n / (i + 1))
        q[k] = round(min(prev, 1.0), 6)
    return q

out = {"tests": ["T1_composition", "T2_position", "T3_matched_window", "T4_recurrent_only"],
       "genes": {}}
for g in GENES:
    gn = GN["genes"].get(g) or {}
    cds = gn.get("cds") or {}
    L = cds.get("protein_aa") or V2["genes"][g].get("protein_length")
    esc = cds.get("nmd_escape_from_codon")
    pn, pf = path[g]["nonsense"], path[g]["frameshift"]
    tol = [v for v in gn.get("plof", []) if v.get("protein_pos") and v.get("lof") == "HC"
           and v.get("pass") and v["consequence"] in ("stop_gained", "frameshift_variant")]
    tn = [v["protein_pos"] for v in tol if v["consequence"] == "stop_gained"]
    tf = [v["protein_pos"] for v in tol if v["consequence"] == "frameshift_variant"]
    rn = [v["protein_pos"] for v in tol if v["consequence"] == "stop_gained" and v["ac"] >= 2]
    rf = [v["protein_pos"] for v in tol if v["consequence"] == "frameshift_variant" and v["ac"] >= 2]
    rec = {"literature_symbol": gn.get("literature_symbol"), "protein_aa": L,
           "nmd_escape_from_codon": esc,
           "n_path_nonsense": len(pn), "n_path_frameshift": len(pf),
           "n_gnomad_nonsense": len(tn), "n_gnomad_frameshift": len(tf),
           "path_fs_share": round(len(pf) / (len(pn) + len(pf)), 3) if (pn or pf) else None,
           "gnomad_fs_share": round(len(tf) / (len(tn) + len(tf)), 3) if (tn or tf) else None}
    if (pn or pf) and (tn or tf):
        rec["T1_p"] = float(fisher_exact([[len(pf), len(pn)], [len(tf), len(tn)]])[1])
    if len(pn) >= 3 and len(pf) >= 3:
        rec["T2_p"] = float(mannwhitneyu(pn, pf, alternative="two-sided")[1])
        rec["T2_median_nonsense_rel"] = round(sorted(pn)[len(pn) // 2] / L, 3) if L else None
        rec["T2_median_frameshift_rel"] = round(sorted(pf)[len(pf) // 2] / L, 3) if L else None
    if esc:
        a = sum(1 for x in pf if x >= esc); b = sum(1 for x in pn if x >= esc)
        c = sum(1 for x in tf if x >= esc); d = sum(1 for x in tn if x >= esc)
        rec.update({"win_path_fs": a, "win_path_ns": b, "win_gnomad_fs": c, "win_gnomad_ns": d})
        if (a + b) and (c + d):
            rec["T3_p"] = float(fisher_exact([[a, b], [c, d]])[1])
    if (pn or pf) and (rn or rf):
        rec["n_gnomad_recurrent_nonsense"] = len(rn)
        rec["n_gnomad_recurrent_frameshift"] = len(rf)
        rec["gnomad_recurrent_fs_share"] = round(len(rf) / (len(rn) + len(rf)), 3)
        rec["T4_p"] = float(fisher_exact([[len(pf), len(pn)], [len(rf), len(rn)]])[1])
    out["genes"][g] = rec

for t in ("T1_p", "T2_p", "T3_p", "T4_p"):
    q = bh([(g, out["genes"][g][t]) for g in GENES if out["genes"][g].get(t) is not None])
    for g, v in q.items():
        out["genes"][g][t.replace("_p", "_q")] = v
json.dump(out, open(os.path.join(HERE, "derived_frame.json"), "w"), indent=1, sort_keys=True)

def nm(g):
    return out["genes"][g].get("literature_symbol") or g
print(f"{'gene':9s} {'ns':>4s} {'fs':>4s} {'fsSh':>5s} | {'gNs':>4s} {'gFs':>4s} {'gFsSh':>6s} "
      f"{'T1q':>8s} | {'T2q':>8s} {'medNs':>6s} {'medFs':>6s} | {'T3q':>8s} | {'T4q':>8s}")
rows = sorted(GENES, key=lambda g: -(out["genes"][g]["n_path_nonsense"] + out["genes"][g]["n_path_frameshift"]))
for g in rows:
    r = out["genes"][g]
    if r["n_path_nonsense"] + r["n_path_frameshift"] == 0:
        continue
    f = lambda k: (f"{r[k]:.1e}" if r.get(k) is not None else "na")
    print(f"{nm(g):9s} {r['n_path_nonsense']:4d} {r['n_path_frameshift']:4d} "
          f"{str(r['path_fs_share']):>5s} | {r['n_gnomad_nonsense']:4d} {r['n_gnomad_frameshift']:4d} "
          f"{str(r['gnomad_fs_share']):>6s} {f('T1_q'):>8s} | {f('T2_q'):>8s} "
          f"{str(r.get('T2_median_nonsense_rel')):>6s} {str(r.get('T2_median_frameshift_rel')):>6s} | "
          f"{f('T3_q'):>8s} | {f('T4_q'):>8s}")
for t in ("T1", "T2", "T3", "T4"):
    s = [nm(g) for g in GENES if (out["genes"][g].get(f"{t}_q") or 1) < 0.05]
    print(f"{t}: {len(s)} genes at q < 0.05: {', '.join(s) if s else 'none'}")
