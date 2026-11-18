import json, os, re, sys, time, urllib.error, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = json.load(open(os.path.join(HERE, "derived_v2.json")))
OUT = os.path.join(HERE, "gnomad_plof.json")
API = "https://gnomad.broadinstitute.org/api"

GENES = sorted(g for g, e in V2["genes"].items() if e["tier"] == "isolated")
DISPLAY = {}
PLOF_CSQ = {"stop_gained", "frameshift_variant", "splice_donor_variant",
            "splice_acceptor_variant", "start_lost", "stop_lost"}
PPOS = re.compile(r'p\.[A-Za-z]{3}(\d+)')

GENE_Q = """{gene(gene_symbol:"%s",reference_genome:GRCh38){
 gene_id symbol chrom
 mane_select_transcript{ensembl_id refseq_id}
 canonical_transcript_id
 variants(dataset:gnomad_r4){variant_id pos consequence hgvsp hgvsc lof lof_filter
  lof_flags flags transcript_id
  exome{ac an homozygote_count filters}
  genome{ac an homozygote_count filters}}}}"""

TX_Q = """{transcript(transcript_id:"%s",reference_genome:GRCh38){
 transcript_id strand chrom start stop exons{feature_type start stop}
 gnomad_constraint{oe_lof oe_lof_upper pLI lof_z exp_lof obs_lof}}}"""


def call(query, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(API, data=json.dumps({"query": query}).encode(),
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=300) as r:
                out = json.load(r)
            if "errors" in out and not out.get("data"):
                raise RuntimeError(json.dumps(out["errors"])[:300])
            return out
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(4 * (i + 1))


def cds_architecture(tx):

    cds = [(e["start"], e["stop"]) for e in tx["exons"] if e["feature_type"] == "CDS"]
    if not cds:
        return None
    cds.sort(reverse=(tx["strand"] == "-"))
    lens = [abs(b - a) + 1 for a, b in cds]
    total = sum(lens)
    last_junction_nt = total - lens[-1]
    penult_50_nt = max(last_junction_nt - 50, 0)
    return {"n_cds_exons": len(cds), "cds_nt": total,
            "protein_aa": total // 3,
            "last_exon_start_codon": last_junction_nt // 3 + 1,
            "nmd_escape_from_codon": penult_50_nt // 3 + 1}


def main():
    out = {"dataset": "gnomad_r4", "genes": {}}
    for i, g in enumerate(GENES, 1):
        sym = DISPLAY.get(g, g)
        try:
            d = call(GENE_Q % sym)["data"]["gene"]
        except Exception as e:
            print(f"[{i}/{len(GENES)}] {sym}: gene query failed: {e}", file=sys.stderr)
            continue
        if not d:
            print(f"[{i}/{len(GENES)}] {sym}: no gene record", file=sys.stderr)
            continue
        mane = (d.get("mane_select_transcript") or {}).get("ensembl_id")
        txid = mane or d.get("canonical_transcript_id")
        arch, constraint = None, None
        if txid:
            try:
                t = call(TX_Q % txid)["data"]["transcript"]
                arch = cds_architecture(t)
                constraint = t.get("gnomad_constraint")
            except Exception as e:
                print(f"  {sym}: transcript query failed: {e}", file=sys.stderr)
        recs = []
        for v in d.get("variants") or []:
            csq = v.get("consequence") or ""
            if csq not in PLOF_CSQ:
                continue
            ex, ge = v.get("exome") or {}, v.get("genome") or {}
            ac = (ex.get("ac") or 0) + (ge.get("ac") or 0)
            an = max(ex.get("an") or 0, ge.get("an") or 0)
            hom = (ex.get("homozygote_count") or 0) + (ge.get("homozygote_count") or 0)
            pass_only = not (ex.get("filters") or ge.get("filters"))
            m = PPOS.search(v.get("hgvsp") or "")
            recs.append({"variant_id": v["variant_id"], "consequence": csq,
                         "hgvsp": v.get("hgvsp"), "lof": v.get("lof"),
                         "lof_filter": v.get("lof_filter"),
                         "flags": v.get("flags") or [], "ac": ac, "an": an,
                         "homozygote_count": hom, "pass": pass_only,
                         "protein_pos": int(m.group(1)) if m else None,
                         "transcript_id": v.get("transcript_id")})
        out["genes"][sym] = {"clinvar_symbol": g, "gene_id": d["gene_id"],
                             "mane_transcript": mane, "transcript_used": txid,
                             "cds": arch, "constraint": constraint,
                             "n_plof": len(recs), "plof": recs}
        hc = sum(1 for r in recs if r["lof"] == "HC")
        print(f"[{i}/{len(GENES)}] {sym}: {len(recs)} pLoF ({hc} LOFTEE HC), "
              f"MANE {mane}, CDS {arch['protein_aa'] if arch else '?'} aa", flush=True)
    json.dump(out, open(OUT, "w"), indent=1, sort_keys=True)
    tot = sum(v["n_plof"] for v in out["genes"].values())
    print(f"\ngenes pulled: {len(out['genes'])}, pLoF variants: {tot}")


if __name__ == "__main__":
    main()
