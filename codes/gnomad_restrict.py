import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lcr_intervals import load, inside

src, dst = sys.argv[1], sys.argv[2]
args = sys.argv[3:]
mane = "--mane" in args
bed = args[args.index("--lcr-bed") + 1] if "--lcr-bed" in args else None
flag = "--lcr-flag" in args
if bed and flag:
    raise SystemExit("give --lcr-bed or --lcr-flag, not both")
if not (mane or bed or flag):
    raise SystemExit("nothing to restrict: give --mane, --lcr-bed BED_GZ, --lcr-flag or a combination")
iv = load(bed) if bed else None
d = json.load(open(src))
kept = dropped = dropped_tx = dropped_lcr = 0
for g, rec in d["genes"].items():
    used = rec.get("transcript_used")
    keep = []
    for v in rec.get("plof", []):
        off_tx = mane and v.get("transcript_id") != used
        ch, pos = v["variant_id"].split("-")[:2]
        in_lcr = (iv is not None and inside(iv, ch, int(pos))) or (flag and "lcr" in (v.get("flags") or []))
        dropped_tx += off_tx
        dropped_lcr += in_lcr
        if not (off_tx or in_lcr):
            keep.append(v)
    dropped += len(rec.get("plof", [])) - len(keep)
    kept += len(keep)
    rec["plof"] = keep
    rec["n_plof"] = len(keep)
d["restricted"] = {"to_transcript_used": mane, "outside_low_complexity_file": os.path.basename(bed) if bed else None,
                   "outside_gnomad_lcr_flag": flag,
                   "records_kept": kept, "records_dropped": dropped,
                   "dropped_off_transcript": dropped_tx, "dropped_in_low_complexity": dropped_lcr}
json.dump(d, open(dst, "w"))
print(json.dumps(d["restricted"]))
