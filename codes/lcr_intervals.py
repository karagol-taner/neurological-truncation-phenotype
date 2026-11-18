import bisect, collections, gzip, json, sys


def load(bed):
    iv = collections.defaultdict(list)
    with gzip.open(bed, "rt") as fh:
        for line in fh:
            c, s, e = line.split()[:3]
            iv[c.replace("chr", "")].append((int(s), int(e)))
    out = {}
    for c, xs in iv.items():
        xs.sort()
        out[c] = ([s for s, _ in xs], [e for _, e in xs])
    return out


def inside(iv, chrom, pos):
    c = str(chrom).replace("chr", "")
    if c not in iv:
        return False
    starts, ends = iv[c]
    i = bisect.bisect_left(starts, pos) - 1
    return i >= 0 and pos <= ends[i]


if __name__ == "__main__":
    bed, gnp, out = sys.argv[1:4]
    iv = load(bed)
    d = json.load(gzip.open(gnp, "rt") if gnp.endswith(".gz") else open(gnp))
    tab = collections.Counter()
    disagree = []
    for sym, rec in d["genes"].items():
        for v in rec.get("plof", []):
            ch, pos = v["variant_id"].split("-")[:2]
            a = inside(iv, ch, int(pos))
            b = "lcr" in (v.get("flags") or [])
            tab[(a, b)] += 1
            if a != b and len(disagree) < 50:
                disagree.append({"gene": sym, "variant_id": v["variant_id"], "in_file": a, "gnomad_flag": b})
    res = {"records": sum(tab.values()),
           "in_file_and_flagged": tab[(True, True)], "in_file_not_flagged": tab[(True, False)],
           "flagged_not_in_file": tab[(False, True)], "neither": tab[(False, False)],
           "first_disagreements": disagree}
    json.dump(res, open(out, "w"), indent=1)
    print({k: v for k, v in res.items() if k != "first_disagreements"})
