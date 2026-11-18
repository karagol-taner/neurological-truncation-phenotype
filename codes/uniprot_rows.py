import importlib.util, json, os, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "uniprot_rows.json")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def query(s):
    q = f"gene:{s} AND organism_id:9606 AND reviewed:true"
    url = ("https://rest.uniprot.org/uniprotkb/search?query="
           + urllib.parse.quote(q)
           + "&fields=accession,length,gene_primary&format=tsv&size=5")
    last = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return {"rows": r.read().decode().splitlines()[1:], "error": None}
        except Exception as e:
            last = repr(e)
            time.sleep(3 * (attempt + 1))
    return {"rows": None, "error": last}


def main(paths):
    symbols = set()
    for i, p in enumerate(paths):
        m = load(p, f"pv2_{i}")
        symbols |= set(m.panel())
    cache = json.load(open(OUT)) if os.path.exists(OUT) else {"symbols": {}}
    todo = sorted(s for s in symbols if s not in cache["symbols"])
    with ThreadPoolExecutor(max_workers=6) as ex:
        for s, res in zip(todo, ex.map(query, todo)):
            if res["error"]:
                json.dump(cache, open(OUT, "w"), indent=0, sort_keys=True)
                raise SystemExit(f"UniProt lookup failed for {s}: {res['error']}")
            cache["symbols"][s] = res
    cache["queried_utc"] = cache.get("queried_utc") or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    json.dump(cache, open(OUT, "w"), indent=0, sort_keys=True)
    print(f"symbols requested {len(symbols)}, newly queried {len(todo)}, cached {len(cache['symbols'])}")


if __name__ == "__main__":
    main(sys.argv[1:])
