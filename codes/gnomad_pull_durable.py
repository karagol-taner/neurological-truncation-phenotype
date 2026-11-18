import gzip, hashlib, importlib.util, json, os, re, socket, sys, time, urllib.error, urllib.request

PULL = os.path.abspath(sys.argv[1])
RAW = os.path.join(PULL, "raw")
STATUS = os.path.join(PULL, "status.json")
NOREC = os.path.join(PULL, "no_gene_record.json")
os.makedirs(RAW, exist_ok=True)
GENE_RE = re.compile(r'gene\(gene_symbol:"([^"]+)"')
TX_RE = re.compile(r'transcript\(transcript_id:"([^"]+)"')
WAITS = (10, 30, 60, 120, 240, 480)
MIN_GAP = 1.0
_last = [0.0]


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True)
    os.replace(tmp, path)


def status(**kw):
    s = json.load(open(STATUS)) if os.path.exists(STATUS) else {}
    s.update(kw)
    s["updated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    write_json(STATUS, s)


def key_of(query):
    m = GENE_RE.search(query)
    if m:
        return f"gene_{m.group(1)}"
    m = TX_RE.search(query)
    if m:
        return f"tx_{m.group(1)}"
    return "q_" + hashlib.sha1(query.encode()).hexdigest()


def fatal(query, err):
    status(state="fatal", failed_query=key_of(query), error=str(err)[:1000])
    print(f"FATAL {key_of(query)}: {err}", file=sys.stderr, flush=True)
    raise SystemExit(2)


def durable_call(query, tries=None):
    k = key_of(query)
    path = os.path.join(RAW, k + ".json.gz")
    if os.path.exists(path):
        with gzip.open(path, "rt") as fh:
            stored = json.load(fh)
        if stored["query"] != query:
            fatal(query, "stored answer was made for a different query")
        return stored["response"]
    last = None
    for attempt in range(len(WAITS) + 1):
        gap = time.time() - _last[0]
        if gap < MIN_GAP:
            time.sleep(MIN_GAP - gap)
        _last[0] = time.time()
        try:
            req = urllib.request.Request(m.API, data=json.dumps({"query": query}).encode(),
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=300) as r:
                out = json.load(r)
            if "errors" in out and not out.get("data"):
                raise RuntimeError(json.dumps(out["errors"])[:300])
            tmp = path + ".tmp"
            with gzip.open(tmp, "wt") as fh:
                json.dump({"query": query, "retrieved_utc": time.strftime(
                    "%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "response": out}, fh)
            os.replace(tmp, path)
            if k.startswith("gene_") and not (out.get("data") or {}).get("gene"):
                nr = json.load(open(NOREC)) if os.path.exists(NOREC) else {}
                nr[k[5:]] = (out.get("errors") or [{}])[0].get("message", "no gene record")
                write_json(NOREC, nr)
            status(state="running", last_query=k,
                   answers_stored=sum(1 for f in os.listdir(RAW) if f.endswith(".json.gz")))
            return out
        except (urllib.error.URLError, urllib.error.HTTPError, socket.timeout, TimeoutError,
                ConnectionError, json.JSONDecodeError, RuntimeError, OSError) as e:
            last = e
            if isinstance(e, urllib.error.HTTPError) and e.code in (400, 404):
                fatal(query, f"HTTP {e.code}")
            if attempt < len(WAITS):
                print(f"retry {attempt + 1} for {k} after {WAITS[attempt]} s: {e}", file=sys.stderr, flush=True)
                time.sleep(WAITS[attempt])
    fatal(query, last)


spec = importlib.util.spec_from_file_location("gnomad_pull", os.path.join(PULL, "gnomad_pull.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
union = json.load(open(os.path.join(PULL, "union_genes.json")))
m.GENES = sorted(union["genes"])
m.call = durable_call
status(state="running", genes_total=len(m.GENES), started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
       failed_query=None, error=None)
m.main()
got = json.load(open(m.OUT))
nr = json.load(open(NOREC)) if os.path.exists(NOREC) else {}
missing = sorted(set(m.GENES) - set(got["genes"]) - set(nr))
if missing:
    status(state="fatal", error=f"genes absent from gnomad_plof.json without a stored reason: {missing[:20]}")
    raise SystemExit(2)
status(state="done", genes_pulled=len(got["genes"]), genes_without_record=sorted(nr),
       genes_without_cds=sorted(g for g, r in got["genes"].items() if not r.get("cds")),
       plof_variants=sum(r["n_plof"] for r in got["genes"].values()))
print("DONE", len(got["genes"]), "genes", flush=True)
