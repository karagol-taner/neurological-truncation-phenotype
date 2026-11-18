import glob
import json
import os
import sys

gate, out = sys.argv[1:3]
cuts = sorted(glob.glob(os.path.join(gate, "*", "*", "codes", "derived_v2.json")))
if len(cuts) != 40:
    raise SystemExit(f"expected 40 gate cuts, found {len(cuts)}")
syms = set()
for p in cuts:
    for e in json.load(open(p))["genes"].values():
        syms |= set(e["clinvar_symbols"])
json.dump({"symbols": sorted(syms), "cuts": [os.path.relpath(p, gate) for p in cuts]}, open(out, "w"), indent=1)
print(len(syms), "symbols from", len(cuts), "cuts")
