import csv, gzip, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lcr_intervals import load, inside

bed, src, dst, rep = sys.argv[1:5]
iv = load(bed)
n_in = n_out = 0
dropped_by_symbol = {}
with gzip.open(src, "rt", encoding="utf-8", errors="replace") as fh, \
        gzip.open(dst, "wt", encoding="utf-8", compresslevel=3, newline="") as out:
    rd = csv.DictReader(fh, delimiter="\t")
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(rd.fieldnames)
    for row in rd:
        n_in += 1
        try:
            pos = int(row.get("PositionVCF") or -1)
        except ValueError:
            pos = -1
        if pos <= 0:
            try:
                pos = int(row.get("Start") or -1)
            except ValueError:
                pos = -1
        if pos > 0 and inside(iv, row.get("Chromosome", ""), pos):
            for s in (row.get("GeneSymbol") or "").split(";"):
                dropped_by_symbol[s] = dropped_by_symbol.get(s, 0) + 1
            continue
        w.writerow([row[k] for k in rd.fieldnames])
        n_out += 1
json.dump({"interval_file": os.path.basename(bed), "rows_in": n_in, "rows_kept": n_out,
           "rows_dropped": n_in - n_out, "dropped_by_symbol": dict(sorted(dropped_by_symbol.items()))},
          open(rep, "w"), indent=1)
print("rows in", n_in, "kept", n_out, "dropped", n_in - n_out)
