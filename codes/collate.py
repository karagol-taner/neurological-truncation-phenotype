import json, os, sys

FLOOR = 10
Q = 0.05
DEFS = ("nmd_50nt_rule", "last_exon", "last_10pct", "last_25pct")


def load(codes, f):
    p = os.path.join(codes, f)
    return json.load(open(p)) if os.path.exists(p) else None


def frame(codes, v2, gn):
    fr = load(codes, "derived_frame.json")
    f2 = load(codes, "derived_fsd2.json")
    fq = load(codes, "derived_fsd_q.json")
    iso = sorted(g for g, e in v2["genes"].items() if e["tier"] == "isolated")
    pool = sorted(g for g in iso if v2["genes"][g]["plp_ptv_positioned"] >= FLOOR)
    genes = fr["genes"]
    n = lambda g: genes[g]["n_path_nonsense"] + genes[g]["n_path_frameshift"]
    tested = {t: sorted(g for g in iso if genes[g].get(f"{t}_p") is not None) for t in ("T1", "T2", "T3", "T4")}
    sig = {t: sorted(g for g in iso if (genes[g].get(f"{t}_q") if genes[g].get(f"{t}_q") is not None else 1) < Q)
           for t in ("T1", "T2", "T3", "T4")}
    no_gnomad = sorted(g for g in pool if g not in gn["genes"])
    no_hc = sorted(g for g in pool if g in gn["genes"] and genes[g]["n_gnomad_nonsense"] + genes[g]["n_gnomad_frameshift"] == 0)
    rows = {}
    for g in sorted(set(sig["T1"]) | set(sig["T4"]) | set(sig["T3"])):
        r = genes[g]
        t2 = f2["genes"].get(g, {})
        rows[g] = {
            "positioned_ptv_panel": v2["genes"][g]["plp_ptv_positioned"], "in_frame_pool": g in pool,
            "path_nonsense": r["n_path_nonsense"], "path_frameshift": r["n_path_frameshift"],
            "path_fs_share": r["path_fs_share"], "gnomad_nonsense": r["n_gnomad_nonsense"],
            "gnomad_frameshift": r["n_gnomad_frameshift"], "gnomad_fs_share": r["gnomad_fs_share"],
            "direction": ("frameshift excess" if (r["path_fs_share"] or 0) > (r["gnomad_fs_share"] or 0)
                          else "nonsense excess"),
            "T1_p": r.get("T1_p"), "T1_q": r.get("T1_q"), "T3_p": r.get("T3_p"), "T3_q": r.get("T3_q"),
            "T4_p": r.get("T4_p"), "T4_q": r.get("T4_q"), "T2_p": r.get("T2_p"), "T2_q": r.get("T2_q"),
            "gnomad_recurrent_fs_share": r.get("gnomad_recurrent_fs_share"),
            "neo_tail": {"n_pathogenic": t2.get("n_pathogenic"), "n_population": t2.get("n_population"),
                         "median_pathogenic": t2.get("median_fsd_pathogenic"),
                         "median_population": t2.get("median_fsd_population"),
                         "max_pathogenic": t2.get("max_fsd_pathogenic"),
                         "mwu_p": t2.get("mwu_p"),
                         "mwu_q": fq["mwu"]["genes"].get(g, {}).get("q"),
                         "long": {k: dict(v, q=fq[f"long_{k}"]["genes"].get(g, {}).get("q"))
                                  for k, v in (t2.get("long") or {}).items()}}}
    return {
        "isolated": len(iso), "frame_pool": len(pool),
        "frame_pool_by_frame_dependence_count": sum(1 for g in iso if n(g) >= FLOOR),
        "tested": {t: len(v) for t, v in tested.items()},
        "frame_pool_tested_T1": sum(1 for g in pool if g in tested["T1"]),
        "frame_pool_failed_on_counts": {"no_gnomad_record": no_gnomad, "no_gnomad_hc_plof": no_hc},
        "isolated_below_floor": len(iso) - len(pool),
        "q_lt_0.05": {t: v for t, v in sig.items()},
        "q_lt_0.05_in_frame_pool": {t: [g for g in v if g in pool] for t, v in sig.items()},
        "separating_genes": rows,
        "neo_tail_q_lt_0.05": {k: fq[k]["q_lt_0.05"] for k in fq},
        "neo_tail_tested": {k: fq[k]["n_genes_tested"] for k in fq},
    }


def pooled(codes):
    ex = load(codes, "derived_fsd_exclusion.json")
    pe = load(codes, "derived_pooled_exclusions.json")
    out = {"neo_tail": {k: v for k, v in ex.items()}}
    nt = pe["neo_tail"]
    out["neo_tail"]["carrying_gene"] = nt.get("carrying_gene")
    if nt.get("carrying_gene"):
        out["neo_tail"][f"pooled_excluding_{nt['carrying_gene']}"] = nt[f"pooled_excluding_{nt['carrying_gene']}"]
    for src in ("windows_submission_level", "windows_aggregated_list"):
        out[src] = {}
        for w, rec in (pe.get(src) or {}).items():
            out[src][w] = {k: v for k, v in rec.items() if k != "leave_one_out_p"}
    return out


def position(codes, suffix):
    ph = load(codes, f"derived_phenotypes{suffix}.json")
    ls = load(codes, f"derived_labels_submission{suffix}.json")
    pq = load(codes, f"derived_position_q{suffix}.json")
    pw = load(codes, f"derived_pilot_windows{suffix}.json")
    pe = load(codes, f"derived_pooled_exclusions{suffix}.json")
    pool = sorted(g for g, r in ph["dual_phenotype_genes"].items() if r["positioned_ptv"] >= FLOOR)
    G = ls["genes"]
    both = sorted(g for g in pool if G[g]["n_area"] and G[g]["n_other"])
    one = sorted(g for g in pool if g not in both)
    res = {"pool": len(pool), "label_buckets": ls["label_buckets"],
           "tested_area_vs_other": len(both),
           "failed_on_counts_one_side_unlabelled": one,
           "tested_area_vs_gnomad": sum(1 for g in pool if G[g]["n_area"] and G[g]["n_gnomad_hc_pass"]),
           "windows": {}, "aggregated_list_q_lt_0.05": {}}
    for w in DEFS:
        sl = pq["submission_level"].get(w, {})
        res["windows"][w] = {c: {"n_genes_tested": sl.get(c, {}).get("n_genes_tested"),
                                 "q_lt_0.05": sl.get(c, {}).get("q_lt_0.05", [])}
                             for c in ("p_area_vs_other", "p_area_vs_gnomad")}
        res["windows"][w]["separating"] = {
            g: dict(G[g]["windows"][w], n_area=G[g]["n_area"], n_other=G[g]["n_other"],
                    q_area_vs_other=sl["p_area_vs_other"]["genes"][g]["q"])
            for g in sl.get("p_area_vs_other", {}).get("q_lt_0.05", [])}
        ag = pq.get("aggregated_list", {}).get(w, {})
        res["aggregated_list_q_lt_0.05"][w] = ag.get("p_area_vs_other", {}).get("q_lt_0.05", [])
    res["pooled"] = {w: {k: v for k, v in rec.items() if k != "leave_one_out_p"}
                     for w, rec in (pe.get("windows_submission_level") or {}).items()}
    res["nominal_p_lt_0.05_area_vs_other"] = sorted(
        g for g in both if any((G[g]["windows"].get(w, {}).get("p_area_vs_other") or 1) < 0.05 for w in DEFS))
    return res


def unit(run):
    codes = os.path.join(run, "codes")
    v2 = load(codes, "derived_v2.json")
    gn = load(codes, "gnomad_plof.json")
    ph = load(codes, "derived_phenotypes.json")
    phs = load(codes, "derived_phenotypes_strictdual.json")
    iso = {g: e for g, e in v2["genes"].items() if e["tier"] == "isolated"}
    pool = {g for g, e in iso.items() if e["plp_ptv_positioned"] >= FLOOR}
    both = {g for g in pool if iso[g]["plp_nonsense"] >= 1 and iso[g]["plp_frameshift"] >= 1}
    dual_d = {g for g, r in ph["dual_phenotype_genes"].items() if r["positioned_ptv"] >= FLOOR}
    dual_s = {g for g, r in phs["dual_phenotype_genes"].items() if r["positioned_ptv"] >= FLOOR}
    ranked = sorted(iso, key=lambda g: (-iso[g]["plp_ptv"], g))
    tot = v2["isolated_ptv_total"]
    run_info = load(codes, "run_panel.json")
    return {
        "run": run_info,
        "gate": {"genes": v2["panel_genes_after_alias_merge"], "isolated": v2["n_isolated"],
                 "syndromic_only": v2["n_syndromic_only"], "rows_grch38": v2["clinvar_rows_matched_grch38"],
                 "ge10_positioned": len(pool), "ge10_with_both_classes": len(both),
                 "ge10_both_and_second_condition_broad_rule": len(both & dual_d),
                 "ge10_both_and_second_condition_strict_rule": len(both & dual_s),
                 "position_pool_broad_rule": len(dual_d), "position_pool_strict_rule": len(dual_s),
                 "largest_gene": ranked[0], "largest_gene_share": round(iso[ranked[0]]["plp_ptv"] / tot, 3),
                 "ATM_in_isolated_tier": "ATM" in iso,
                 "ATM_share": round(iso["ATM"]["plp_ptv"] / tot, 3) if "ATM" in iso else None},
        "frame": frame(codes, v2, gn),
        "pooled": pooled(codes),
        "position_broad_rule": position(codes, ""),
        "position_strict_rule": position(codes, "_strictdual"),
    }


out = {}
for arg in sys.argv[2:]:
    name, run = arg.split("=", 1)
    out[name] = unit(run)
json.dump(out, open(sys.argv[1], "w"), indent=1, sort_keys=True)
for name, u in out.items():
    g, f = u["gate"], u["frame"]
    print(f"\n#### {name}: genes {g['genes']} isolated {g['isolated']} | >=10 {g['ge10_positioned']} both {g['ge10_with_both_classes']} "
          f"dual(broad) {g['ge10_both_and_second_condition_broad_rule']} dual(strict) {g['ge10_both_and_second_condition_strict_rule']} "
          f"| largest {g['largest_gene']} {g['largest_gene_share']}")
    print(f"  frame: T1 tested {f['tested']['T1']} (pool {f['frame_pool_tested_T1']}/{f['frame_pool']}); q<0.05 "
          + "; ".join(f"{t} {len(v)}" for t, v in f["q_lt_0.05"].items()))
    for gname, r in f["separating_genes"].items():
        nt = r["neo_tail"]
        print(f"   {gname:9s} pool={r['in_frame_pool']!s:5s} ns {r['path_nonsense']:4d} fs {r['path_frameshift']:4d} share {r['path_fs_share']} "
              f"vs gnomAD {r['gnomad_fs_share']} ({r['gnomad_nonsense']}/{r['gnomad_frameshift']}) T1 q {r['T1_q']} T4 q {r['T4_q']} T3 q {r['T3_q']} "
              f"| tail med {nt['median_pathogenic']} vs {nt['median_population']} q {nt['mwu_q']}")
    nt = u["pooled"]["neo_tail"]
    for k in sorted(k for k in nt if k.startswith("pooled")):
        print(f"  neo-tail {k}: n {nt[k]['n_path']}/{nt[k]['n_pop']} med {nt[k]['med_path']}/{nt[k]['med_pop']} p {nt[k].get('mwu_p')}")
    for rule in ("position_broad_rule", "position_strict_rule"):
        p = u[rule]
        print(f"  {rule}: pool {p['pool']}, tested {p['tested_area_vs_other']}, one side unlabelled {len(p['failed_on_counts_one_side_unlabelled'])}; "
              + "; ".join(f"{w} q<0.05 {p['windows'][w]['p_area_vs_other']['q_lt_0.05']}" for w in DEFS))
        for w in DEFS:
            pl = p["pooled"].get(w, {})
            print("     pooled " + w + ": " + ", ".join(f"{k} p={v['p']:.2e}" for k, v in pl.items() if k.startswith("pooled"))
                  + f" (carrying {pl.get('carrying_gene')})")
