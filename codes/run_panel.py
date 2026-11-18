import argparse, collections, csv, importlib.util, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = json.load(open(os.path.join(HERE, "uniprot_rows.json")))["symbols"]


def parse_first(rows):
    if not rows:
        return None
    f = rows[0].split("\t")
    if len(f) >= 3:
        try:
            return {"acc": f[0], "length": int(f[1]), "primary": f[2].strip()}
        except ValueError:
            return None
    return None


def parse_strict(s, rows):
    hits = []
    for r in rows or []:
        f = r.split("\t")
        if len(f) >= 3 and s in [x.strip() for x in f[2].split(";")]:
            try:
                hits.append({"acc": f[0], "length": int(f[1]), "primary": f[2].strip()})
            except ValueError:
                pass
    return hits


def grouping_names(gcs_path, n):
    genes = collections.defaultdict(set)
    with open(gcs_path, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = (r.get("AssociatedGenes") or "").strip()
            if s:
                genes[(r.get("DiseaseName") or "").strip()].add(s)
    return {d: len(g) for d, g in genes.items() if len(g) >= n}


def filtered_map(gcs_path, n, out_path, scope="names", in_area=None):

    drop = grouping_names(gcs_path, n)
    lines = open(gcs_path, encoding="utf-8", errors="replace").read().splitlines(keepends=True)
    header, body = lines[0], lines[1:]
    cols = header.lstrip("#").rstrip("\n").split("\t")
    i, j = cols.index("DiseaseName"), cols.index("AssociatedGenes")
    field = lambda f, k: f[k].strip() if len(f) > k else ""
    routes_only = set()
    if scope == "routes":
        names = collections.defaultdict(set)
        for line in body:
            f = line.rstrip("\n").split("\t")
            if field(f, j) and in_area(field(f, i)):
                names[field(f, j)].add(field(f, i))
        routes_only = {s for s, ns in names.items() if ns and all(x in drop for x in ns)}
    kept = removed = 0
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(header)
        for line in body:
            f = line.rstrip("\n").split("\t")
            d, s = field(f, i), field(f, j)
            if (scope == "names" and d in drop) or \
                    (scope == "routes" and s in routes_only and in_area(d)):
                removed += 1
                continue
            out.write(line)
            kept += 1
    return drop, kept, removed, sorted(routes_only)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("codes_dir")
    ap.add_argument("--resolution", choices=("first", "strict"), required=True)
    ap.add_argument("--grouping", type=int, default=0)
    ap.add_argument("--grouping-scope", choices=("names", "routes"), default="names")
    a = ap.parse_args()
    codes = os.path.abspath(a.codes_dir)
    spec = importlib.util.spec_from_file_location("panel_v2", os.path.join(codes, "panel_v2.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    info = {"resolution": a.resolution, "grouping_threshold": a.grouping,
            "gene_condition_source_id": os.path.realpath(m.GCS),
            "variant_summary": os.path.realpath(m.VS)}
    if a.grouping:
        fpath = os.path.join(codes, f"gene_condition_source_id.grouping_{a.grouping_scope}_lt{a.grouping}")
        in_area = lambda d: bool(m.DISEASE.search(d)) and not m.EXCL.search(d)
        drop, kept, removed, routes_only = filtered_map(m.GCS, a.grouping, fpath,
                                                        a.grouping_scope, in_area)
        m.GCS = fpath
        info.update({"grouping_scope": a.grouping_scope, "grouping_names": len(drop),
                     "grouping_names_in_area": sorted(d for d in drop if in_area(d)),
                     "symbols_whose_only_route_is_a_grouping_name": routes_only,
                     "map_rows_kept": kept, "map_rows_removed": removed})
    info["map_used"] = os.path.abspath(m.GCS)
    failures = []

    def resolve(symbols):
        out = {}
        for s in symbols:
            c = ROWS.get(s)
            if c is None or c["error"]:
                raise SystemExit(f"no cached UniProt answer for {s}")
            if a.resolution == "first":
                r = parse_first(c["rows"])
                if r:
                    out[s] = r
            else:
                hits = parse_strict(s, c["rows"])
                if len(hits) == 1:
                    out[s] = hits[0]
                else:
                    failures.append([s, len(hits)])
        return out

    m.uniprot = resolve
    m.main()
    info["resolution_failures"] = sorted(failures)
    json.dump(sorted(failures), open(os.path.join(codes, "resolution_failures.json"), "w"))
    json.dump(info, open(os.path.join(codes, "run_panel.json"), "w"), indent=1, sort_keys=True)


if __name__ == "__main__":
    main()
