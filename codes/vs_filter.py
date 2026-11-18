import csv, gzip, json, sys
syms = set(json.load(open(sys.argv[1]))["symbols"])
n_in = n_out = 0
with gzip.open(sys.argv[2], "rt", encoding="utf-8", errors="replace") as fh, \
        gzip.open(sys.argv[3], "wt", encoding="utf-8", compresslevel=3, newline="") as out:
    rd = csv.DictReader(fh, delimiter="\t")
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(rd.fieldnames)
    for row in rd:
        n_in += 1
        if row.get("Assembly") == "GRCh38" and syms & set((row.get("GeneSymbol") or "").split(";")):
            w.writerow([row[k] for k in rd.fieldnames]); n_out += 1
print("rows in", n_in, "rows kept", n_out)
