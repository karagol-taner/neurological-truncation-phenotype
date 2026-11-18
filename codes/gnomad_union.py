import json, os, sys

GATE, PULL = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
UNITS = ("brain_whole", "brain_epilepsy", "brain_neurodevelopmental",
         "brain_neurodegeneration", "brain_malformation")


def v2(variant, unit):
    return json.load(open(os.path.join(GATE, variant, unit, "codes", "derived_v2.json")))


src, steps = {}, []
for variant in ("screen", "strict"):
    for unit in UNITS:
        d = v2(variant, unit)
        iso = sorted(g for g, e in d["genes"].items() if e["tier"] == "isolated")
        new = [g for g in iso if g not in src]
        for g in iso:
            src.setdefault(g, []).append(f"{variant}/{unit}/isolated")
        steps.append({"step": f"{variant}/{unit} isolated tier", "isolated": len(iso),
                      "added": len(new), "added_genes": new if len(new) <= 50 else None})
unmerged = set()
for unit in UNITS:
    first = v2("screen", unit)
    strict = v2("strict", unit)
    merged_first = {s for k, v in first["alias_merges"].items() for s in v}
    merged_strict = {s for k, v in strict["alias_merges"].items() for s in v}
    unmerged |= {s for s in merged_first - merged_strict if s in strict["genes"]}
new = sorted(g for g in unmerged if g not in src)
for g in sorted(unmerged):
    src.setdefault(g, []).append("strict/unmerged")
steps.append({"step": "symbols the strict resolution unmerges", "unmerged": sorted(unmerged),
              "added": len(new), "added_genes": new})
os.makedirs(PULL, exist_ok=True)
json.dump({"n_genes": len(src), "steps": steps, "genes": src},
          open(os.path.join(PULL, "union_genes.json"), "w"), indent=1, sort_keys=True)
for s in steps:
    print(s["step"], "| added", s["added"], s.get("added_genes") or "")
print("union:", len(src))
