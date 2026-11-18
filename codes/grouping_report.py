import collections, csv, importlib.util, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
GATE, OUT = os.path.abspath(sys.argv[1]), sys.argv[2]
GCS = os.path.join(os.path.dirname(HERE), "gene_condition_source_id")
UNITS = ("brain_whole", "brain_epilepsy", "brain_neurodevelopmental", "brain_neurodegeneration", "brain_malformation")


def v2(variant, unit):
    return json.load(open(os.path.join(GATE, variant, unit, "codes", "derived_v2.json")))


def frame(d):
    return {g for g, e in d["genes"].items() if e["tier"] == "isolated" and e["plp_ptv_positioned"] >= 10}


name_genes = collections.defaultdict(set)
for r in csv.DictReader(open(GCS, encoding="utf-8", errors="replace"), delimiter="\t"):
    s = (r.get("AssociatedGenes") or "").strip()
    if s:
        name_genes[(r.get("DiseaseName") or "").strip()].add(s)
spec = importlib.util.spec_from_file_location("pv", os.path.join(GATE, "screen", "brain_whole", "codes", "panel_v2.py"))
pv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)
in_area = lambda c: bool(pv.DISEASE.search(c)) and not pv.EXCL.search(c)
out = {"names": {}, "only_route_screen_panel": {}, "frame_lost_whole": {}, "re_tiered_whole": {}}
for thr in (10, 20):
    out["names"][f"ge{thr}"] = {"in_map": sum(1 for g in name_genes.values() if len(g) >= thr),
                                "in_area": {d: len(g) for d, g in name_genes.items() if len(g) >= thr and in_area(d)}}
    d = v2("screen", "brain_whole")
    fr = frame(d)
    only = sorted(g for g, e in d["genes"].items()
                  if (set(e["isolated"]) | set(e["syndromic"]))
                  and all(len(name_genes[n]) >= thr for n in set(e["isolated"]) | set(e["syndromic"])))
    out["only_route_screen_panel"][f"ge{thr}"] = {
        "genes": len(only), "isolated": sum(1 for g in only if d["genes"][g]["tier"] == "isolated"),
        "frame_pool": sorted(g for g in only if g in fr), "all": only}
    base = v2("strict_dicer", "brain_whole")
    for scope, var in (("routes", f"strict_dicer_r{thr}"), ("names", f"strict_dicer_g{thr}")):
        new = v2(var, "brain_whole")
        out["frame_lost_whole"][f"{scope}_ge{thr}"] = sorted(frame(base) - frame(new))
        out["re_tiered_whole"][f"{scope}_ge{thr}"] = sorted(
            g for g in set(base["genes"]) & set(new["genes"]) if base["genes"][g]["tier"] != new["genes"][g]["tier"])
hpe = "Holoprosencephaly sequence"
d = v2("screen", "brain_whole")
out["holoprosencephaly_sequence"] = {"genes_in_map": len(name_genes[hpe]),
                                     "only_route_for": sorted(g for g, e in d["genes"].items()
                                                              if set(e["isolated"]) | set(e["syndromic"]) == {hpe})}
out["genetic_hearing_loss_genes_in_map"] = len(name_genes["Genetic hearing loss"])
nd = v2("screen", "brain_neurodegeneration")["genes"].get("ATM", {})
out["ATM_names_in_neurodegeneration_unit"] = sorted(set(nd.get("isolated", [])) | set(nd.get("syndromic", [])))
json.dump(out, open(OUT, "w"), indent=1, sort_keys=True)
print(json.dumps({k: v for k, v in out.items() if k != "only_route_screen_panel"}, indent=1)[:2500])
for k, v in out["only_route_screen_panel"].items():
    print(k, v["genes"], v["isolated"], v["frame_pool"])
