import csv, gzip, sys
ids = set()
with gzip.open(sys.argv[1], "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        ids.add(row["VariationID"])
n_in = n_out = 0
with gzip.open(sys.argv[2], "rt", encoding="utf-8", errors="replace") as fh, \
        gzip.open(sys.argv[3], "wt", encoding="utf-8", compresslevel=3) as out:
    for line in fh:
        if line.startswith("#"):
            out.write(line); continue
        n_in += 1
        if line.split("\t", 1)[0] in ids:
            out.write(line); n_out += 1
print("variation ids", len(ids), "records in", n_in, "kept", n_out)
