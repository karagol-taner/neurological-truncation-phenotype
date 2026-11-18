import collections, gzip, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lcr_intervals import load, inside

bed, gnp, v2p, out = sys.argv[1:5]
iv = load(bed)
gn = json.load(gzip.open(gnp, "rt") if gnp.endswith(".gz") else open(gnp))["genes"]
v2 = json.load(open(v2p))["genes"]
iso = {g for g, e in v2.items() if e["tier"] == "isolated"}
FSTER = re.compile(r"fsTer(\d+)", re.I)

agree = collections.Counter()
for sym, rec in gn.items():
    for v in rec.get("plof", []):
        ch, pos = v["variant_id"].split("-")[:2]
        agree[(inside(iv, ch, int(pos)), "lcr" in (v.get("flags") or []))] += 1

genes, tot = {}, collections.Counter()
for sym, rec in sorted(gn.items()):
    g = rec.get("clinvar_symbol") or sym
    if g not in iso or not rec.get("cds"):
        continue
    c = collections.Counter()
    esc = rec["cds"]["nmd_escape_from_codon"]
    for v in rec["plof"]:
        if not (v.get("lof") == "HC" and v.get("pass") and v.get("protein_pos")
                and v["consequence"] in ("stop_gained", "frameshift_variant")):
            continue
        k = "fs" if v["consequence"] == "frameshift_variant" else "ns"
        ch, pos = v["variant_id"].split("-")[:2]
        lcr = inside(iv, ch, int(pos))
        flag = "lcr" in (v.get("flags") or [])
        c[k] += 1
        c[k + "_lcr"] += lcr
        c[k + "_flag"] += flag
        if v["protein_pos"] >= esc:
            c["win"] += 1
            c["win_lcr"] += lcr
        if k == "fs" and FSTER.search(v.get("hgvsp") or ""):
            c["tails"] += 1
            c["tails_lcr"] += lcr
            c["tails_flag"] += flag
        if k == "fs" and v["ac"] >= 2:
            c["fs_rec"] += 1
            c["fs_rec_lcr"] += lcr
    genes[g] = dict(c)
    tot.update(c)

json.dump({"interval_file": os.path.basename(bed),
           "agreement_all_records": {"records": sum(agree.values()),
                                     "in_file_and_flagged": agree[(True, True)],
                                     "in_file_not_flagged": agree[(True, False)],
                                     "flagged_not_in_file": agree[(False, True)],
                                     "neither": agree[(False, False)]},
           "isolated_tier_totals": dict(tot), "genes": genes}, open(out, "w"), indent=1, sort_keys=True)
print(json.dumps(dict(tot)))
