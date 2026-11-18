#!/usr/bin/env python3

import csv
import json
import math
import os
import re
from collections import defaultdict

import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

from figstyle import plt, tidy, save, RED, S3, GREY, INK, INK2, MUTED

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
FIGS = os.path.join(HERE, "..", "figures")
FLOOR = 10
DEFS = ("nmd_50nt_rule", "last_exon", "last_10pct", "last_25pct")
DEF_LABEL = {"nmd_50nt_rule": "Fifty-nucleotide rule", "last_exon": "Final exon",
             "last_10pct": "Last tenth", "last_25pct": "Last quarter"}


def load(*p):
    return json.load(open(os.path.join(DATA, *p)))


RS = load("results_summary.json")
W = RS["brain_whole"]
V2 = load("derived_v2.json")
FR = load("derived_frame.json")["genes"]
F2 = load("derived_fsd2.json")["genes"]
LS = load("derived_labels_submission.json")["genes"]
PQ = load("derived_position_q.json")["submission_level"]
VT = list(csv.DictReader(open(os.path.join(DATA, "variant_table.csv"))))

ISO = {g: e for g, e in V2["genes"].items() if e["tier"] == "isolated"}
SIG = W["frame"]["q_lt_0.05"]
SEP = sorted(set(SIG["T1"]) | set(SIG["T4"]))

FRL = load("sensitivity", "lcrsym_whole", "derived_frame.json")["genes"]
ROBUST = set(SEP)
for run in ("mane_whole", "lcrsym_whole", "manelcrsym_whole"):
    s_ = RS[run]["frame"]["q_lt_0.05"]
    ROBUST &= set(s_["T1"]) | set(s_["T4"])
SEPG = W["frame"]["separating_genes"]
NS_DIR = [g for g in SEP if SEPG[g]["direction"] == "nonsense excess"]
FS_DIR = [g for g in SEP if SEPG[g]["direction"] == "frameshift excess"]


def label_point(ax, x, y, text, dx=4, dy=3, ha="left", color=INK2, size=7.5):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points",
                ha=ha, va="bottom", fontsize=size, color=color)


def ntot(e, side):
    return (e["n_path_nonsense"] + e["n_path_frameshift"] if side == "p"
            else e["n_gnomad_nonsense"] + e["n_gnomad_frameshift"])


def nwin(e, side):
    return (e["win_path_fs"] + e["win_path_ns"] if side == "p"
            else e["win_gnomad_fs"] + e["win_gnomad_ns"])


def figure1():
    fig = plt.figure(figsize=(7.4, 4.9))
    gsp = fig.add_gridspec(2, 2, width_ratios=[1.1, 1.0], hspace=0.95, wspace=0.75)

    ax = fig.add_subplot(gsp[:, 0])
    ranked = sorted(ISO.items(), key=lambda kv: (-kv[1]["plp_ptv"], kv[0]))
    ranked = [(g, e) for g, e in ranked if e["plp_ptv"] > 0]
    ys = [e["plp_ptv"] for _, e in ranked]
    xs = np.arange(1, len(ys) + 1)
    ax.plot(xs, ys, color=RED, lw=1.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(FLOOR, color=MUTED, lw=0.8, ls="--")
    ax.text(1.2, FLOOR * 0.62, f"{FLOOR} variants", ha="left", va="top", fontsize=7.5, color=INK2)
    total = V2["isolated_ptv_total"]
    g, e = ranked[0]
    ax.plot([1], [e["plp_ptv"]], "o", color=INK, ms=4)
    label_point(ax, 1, e["plp_ptv"], f"{g}, {100.0 * e['plp_ptv'] / total:.1f}% of variants",
                dx=7, dy=-4, size=7.5)
    tidy(ax, ylab="Pathogenic truncating variants", xlab="Gene rank, isolated tier",
         title="a  Variants per gene", grid="both")

    ax = fig.add_subplot(gsp[0, 1])
    ga, fa = W["gate"], W["frame"]
    pd_ = W["position_broad_rule"]
    steps = [("Panel", ga["genes"], GREY),
             ("Isolated tier", ga["isolated"], GREY),
             (f"{FLOOR} or more variants", ga["ge10_positioned"], RED),
             ("Frame pool tested", fa["frame_pool_tested_T1"], RED),
             ("Second condition", ga["position_pool_broad_rule"], INK2),
             ("Position axis tested", pd_["tested_area_vs_other"], INK2)]
    for i, (lab, v, c) in enumerate(steps):
        ax.barh(i, v, color=c, height=0.66)
        ax.text(v + ga["genes"] * 0.02, i, f"{v:,}", va="center", ha="left", fontsize=7.5, color=INK2)
    ax.set_yticks(range(len(steps)))
    ax.set_yticklabels([s[0] for s in steps])
    ax.invert_yaxis()
    ax.set_xlim(0, ga["genes"] * 1.3)
    tidy(ax, xlab="Genes", title="b  Genes reaching each step", grid="x")

    ax = fig.add_subplot(gsp[1, 1])
    buckets = [("area", "Neurological", RED), ("other", "Other", INK2),
               ("both", "Both", MUTED), ("unlabelled", "Neither", GREY)]
    for i, rule in enumerate(("position_broad_rule", "position_strict_rule")):
        lb = W[rule]["label_buckets"]
        tot = sum(lb.values())
        left = 0.0
        for key, lab, col in buckets:
            w = 100.0 * lb[key] / tot
            ax.barh(i, w, left=left, color=col, height=0.6)
            left += w
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Broad rule", "Strict rule"])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    handles = [Rectangle((0, 0), 1, 1, color=c) for _, _, c in buckets]
    ax.legend(handles, [b[1] for b in buckets], ncol=4, loc="upper center",
              bbox_to_anchor=(0.42, -0.42), handlelength=1.0, columnspacing=0.8)
    tidy(ax, xlab="Per cent of variants", title="c  Reported phenotype", grid="x")
    save(fig, "figure1")


def figure2():
    fig = plt.figure(figsize=(7.4, 4.9))
    gsp = fig.add_gridspec(2, 2, hspace=0.8, wspace=0.45)

    ax = fig.add_subplot(gsp[0, 0])
    cols = {"nmd_50nt_rule": RED, "last_exon": INK, "last_10pct": MUTED, "last_25pct": GREY}
    top = 0.0
    for w in DEFS:
        ps = sorted(r["p"] for r in PQ[w]["p_area_vs_other"]["genes"].values()
                    if r.get("p") is not None)
        m = len(ps)
        exp = -np.log10((np.arange(1, m + 1) - 0.5) / m)
        obs = -np.log10(np.array(ps))
        ax.plot(exp, obs, "o", ms=2.6, color=cols[w], label=DEF_LABEL[w], alpha=0.9)
        top = max(top, exp.max(), obs.max())
    ax.plot([0, top * 1.05], [0, top * 1.05], color=MUTED, lw=0.8, ls="--")
    ax.set_xlim(0, top * 1.05)
    ax.set_ylim(0, top * 1.05)
    ax.legend(loc="upper left", handletextpad=0.2, borderaxespad=0.2, fontsize=6.5)
    tidy(ax, xlab="Expected −log₁₀ P", ylab="Observed −log₁₀ P",
         title="a  Neurological against other", grid="both")

    ax = fig.add_subplot(gsp[0, 1])
    pts = []
    for g, r in LS.items():
        w = r["windows"].get("nmd_50nt_rule")
        if not (r["n_area"] and r["n_other"]) or not w:
            continue
        pts.append((g, w["threshold_rel"], r["n_area"] / r["n_other"],
                    r["n_area"] + r["n_other"]))
    carry = W["position_broad_rule"]["pooled"]["nmd_50nt_rule"]["carrying_gene"]
    ax.scatter([p[1] for p in pts], [p[2] for p in pts], s=[6 + 0.08 * p[3] for p in pts],
               color=GREY, edgecolor="none", alpha=0.85)
    for g, x, y, nn in pts:
        if g in (carry, "ATM"):
            ax.scatter([x], [y], s=6 + 0.08 * nn, color=INK2 if g == carry else RED, edgecolor="none")
            label_point(ax, x, y, g, dx=7, dy=-2 if g == carry else -11, color=INK)
    ax.set_yscale("log")
    ax.axhline(1, color=MUTED, lw=0.8, ls="--")
    ax.set_xlim(0, 1.02)
    tidy(ax, xlab="Window opens, fraction of protein", ylab="Neurological to other labels",
         title="b  Strata of a pooled table", grid="both")

    ax = fig.add_subplot(gsp[1, :])
    pooled = W["position_broad_rule"]["pooled"]
    rows_ = []
    for w in DEFS:
        rec = pooled[w]
        c = rec["carrying_gene"]
        rows_.append((DEF_LABEL[w], rec["pooled_all"]["p"], rec["pooled_excluding_ATM"]["p"],
                      rec[f"pooled_excluding_{c}"]["p"], c))
    ymax = max(-math.log10(r[1]) for r in rows_)
    for i, (lab, pa, pat, pc, c) in enumerate(rows_):
        for j, (p, col, mk) in enumerate(((pa, INK, "o"), (pat, RED, "s"), (pc, INK2, "D"))):
            ax.plot([i + (j - 1) * 0.2], [-math.log10(p)], mk, color=col, ms=5)
        ax.annotate(c, (i + 0.2, -math.log10(pc)), xytext=(5, 7), textcoords="offset points",
                    ha="left", va="bottom", fontsize=7, color=INK2)
    ax.axhline(-math.log10(0.05), color=MUTED, lw=0.8, ls="--")
    ax.set_xticks(range(len(rows_)))
    ax.set_xticklabels([r[0] for r in rows_])
    ax.set_xlim(-0.5, len(rows_) - 0.3)
    ax.set_ylim(-ymax * 0.03, ymax * 1.15)
    hs = [Line2D([0], [0], marker=m, color=c, ls="none", ms=5)
          for m, c in (("o", INK), ("s", RED), ("D", INK2))]
    ax.legend(hs, ["All genes", "Without ATM", "Without the carrying gene"], ncol=3,
              loc="upper right", fontsize=7)
    tidy(ax, ylab="−log₁₀ pooled P", title="c  Pooled window test", grid="y")
    save(fig, "figure2")


def figure3():
    genes = W["position_broad_rule"]["windows"]["nmd_50nt_rule"]["p_area_vs_gnomad"]["q_lt_0.05"]
    share = lambda g, k, n: LS[g]["windows"]["nmd_50nt_rule"][k] / LS[g][n]
    pop_more = [g for g in genes if share(g, "gnomad_in", "n_gnomad_hc_pass") > share(g, "area_in", "n_area")]
    order = sorted(pop_more) + sorted(g for g in genes if g not in pop_more)
    pos = defaultdict(lambda: {"b": [], "g": []})
    for r in VT:
        g = r["gene"]
        if g not in order or not r["rel_pos"]:
            continue
        if r["source"] == "clinvar_pathogenic" and r["phenotype_side"] == "area":
            pos[g]["b"].append(float(r["rel_pos"]))
        elif r["source"] == "gnomad_plof":
            pos[g]["g"].append(float(r["rel_pos"]))
    for g in order:
        assert len(pos[g]["b"]) == LS[g]["n_area"], g
        assert len(pos[g]["g"]) == LS[g]["n_gnomad_hc_pass"], g
    qd = PQ["nmd_50nt_rule"]["p_area_vs_gnomad"]["genes"]
    xmax = max(1.0, max(max(v["b"] + v["g"]) for v in pos.values()))
    x_cnt, x_q = xmax + 0.04, xmax + 0.16
    fig, axs = plt.subplots(len(order), 1, figsize=(7.2, 4.8), sharex=True)
    rng = np.random.default_rng(7)
    for ax, g in zip(axs, order):
        w = LS[g]["windows"]["nmd_50nt_rule"]
        thr = w["threshold_rel"]
        ax.add_patch(Rectangle((thr, -0.55), 1.0 - thr, 2.1, color="#e9e8e4", lw=0))
        for row, key, col in ((1, "b", S3), (0, "g", MUTED)):
            xs = np.array(pos[g][key])
            ax.scatter(xs, row + rng.uniform(-0.22, 0.22, len(xs)), s=5, color=col,
                       edgecolor="none", alpha=0.85, clip_on=False)
        ax.text(x_cnt, 1, f"{w['area_in']}/{LS[g]['n_area']}", va="center", ha="left",
                fontsize=7, color=S3, transform=ax.transData)
        ax.text(x_cnt, 0, f"{w['gnomad_in']}/{LS[g]['n_gnomad_hc_pass']}", va="center", ha="left",
                fontsize=7, color=INK2, transform=ax.transData)
        q = qd[g]["q"]
        ax.text(x_q, 0.5, ("q < 0.001" if q < 0.001 else f"q = {q:.3f}"), va="center",
                ha="left", fontsize=7, color=INK, transform=ax.transData)
        ax.set_ylim(-0.6, 1.6)
        ax.set_yticks([0.5])
        ax.set_yticklabels([g])
        ax.tick_params(axis="y", length=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.set_xlim(0, xmax + 0.005)
        ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axs[-1].set_xlabel("Codon position, fraction of protein")
    hs = [Line2D([0], [0], marker="o", color=S3, ls="none", ms=4),
          Line2D([0], [0], marker="o", color=MUTED, ls="none", ms=4),
          Rectangle((0, 0), 1, 1, color="#e9e8e4")]
    fig.legend(hs, ["Pathogenic, neurological condition", "gnomAD", "Escape window"],
               ncol=3, loc="lower center", bbox_to_anchor=(0.45, 0.9), fontsize=7.5,
               handletextpad=0.3, columnspacing=1.2)
    fig.subplots_adjust(hspace=0.3, right=0.78, top=0.88)
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(os.path.join(FIGS, "figure3.png"), dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(FIGS, "figure3.svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote figure3.png and .svg")


def figure4():
    fig = plt.figure(figsize=(7.4, 4.9))
    gsp = fig.add_gridspec(2, 2, width_ratios=[1.1, 1.0], hspace=0.75, wspace=0.5)

    ax = fig.add_subplot(gsp[:, 0])
    pool = sorted(g for g, e in ISO.items() if e["plp_ptv_positioned"] >= FLOOR and FR[g].get("T1_p") is not None)
    for g in pool:
        e = FR[g]
        col = RED if g in NS_DIR else INK if g in FS_DIR else GREY
        ax.scatter([e["gnomad_fs_share"]], [e["path_fs_share"]], s=4 + 2.2 * math.sqrt(ntot(e, "p")),
                   color=col, edgecolor="none", alpha=0.95 if g in SEP else 0.5,
                   zorder=3 if g in SEP else 2)
    offs = {"ARID1A": (7, -6), "ATM": (7, -3), "CIC": (-6, 3), "IQSEC2": (-7, -2),
            "IRF2BPL": (-7, -12), "KMT2B": (6, 2), "SETD1B": (-7, 3), "SRRM2": (7, -6),
            "ZFHX3": (7, -9)}
    offs.update({"CIC": (-18, 14), "SETD1B": (-26, 4), "IRF2BPL": (14, -12), "IQSEC2": (-22, -10)})
    for g in SEP:
        e = FR[g]
        dx, dy = offs.get(g, (5, 3))
        ax.annotate(g, (e["gnomad_fs_share"], e["path_fs_share"]), xytext=(dx, dy),
                    textcoords="offset points", ha="left" if dx > 0 else "right",
                    va="center", fontsize=7, color=INK,
                    arrowprops=(dict(arrowstyle="-", lw=0.5, color=MUTED, shrinkA=0, shrinkB=2)
                                if abs(dx) > 12 else None))
    ax.plot([0, 1], [0, 1], color=MUTED, lw=0.8, ls="--")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.04)
    tidy(ax, xlab="gnomAD frameshift share", ylab="Pathogenic frameshift share",
         title="a  Composition, frame pool", grid="both")

    ax = fig.add_subplot(gsp[0, 1])
    order = sorted(SEP, key=lambda g: (-nwin(FR[g], "g") / ntot(FR[g], "g"), g))
    for i, g in enumerate(order):
        e = FR[g]
        pw = 100.0 * nwin(e, "p") / ntot(e, "p")
        gw = 100.0 * nwin(e, "g") / ntot(e, "g")
        ax.plot([pw, gw], [i, i], color=GREY, lw=1.2, zorder=1)
        ax.plot([gw], [i], "o", color=MUTED, ms=4.5, zorder=2)
        ax.plot([pw], [i], "o", color=RED, ms=4.5, zorder=3)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(-4, 104)
    tidy(ax, xlab="Per cent in escape window", title="b  Escape window", grid="x")
    for lab in ax.get_yticklabels():
        lab.set_fontweight("bold" if lab.get_text() in ROBUST else "normal")

    ax = fig.add_subplot(gsp[1, 1])
    for i, g in enumerate(order):
        e, el = FR[g], FRL[g]
        a, r_, p = e["gnomad_fs_share"], e["gnomad_recurrent_fs_share"], e["path_fs_share"]
        al, pl = el["gnomad_fs_share"], el["path_fs_share"]
        lo, hi = min(a, r_, p, al, pl), max(a, r_, p, al, pl)
        ax.plot([lo, hi], [i, i], color=GREY, lw=1.0, zorder=1)
        ax.plot([al], [i], "o", mfc="white", mec=MUTED, mew=1.0, ms=4.5, zorder=1.5)
        ax.plot([pl], [i], "o", mfc="white", mec=RED, mew=1.0, ms=4.5, zorder=1.6)
        ax.plot([a], [i], "o", color=MUTED, ms=4.5, zorder=2)
        ax.plot([r_], [i], "|", color=INK, ms=8, mew=1.4, zorder=2)
        ax.plot([p], [i], "o", color=RED, ms=4.5, zorder=3)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7)
    for lab in ax.get_yticklabels():
        lab.set_fontweight("bold" if lab.get_text() in ROBUST else "normal")
    ax.invert_yaxis()
    ax.set_xlim(0.35, 1.03)
    hs = [Line2D([0], [0], marker="o", color=RED, ls="none", ms=4.5),
          Line2D([0], [0], marker="o", color=MUTED, ls="none", ms=4.5),
          Line2D([0], [0], marker="|", color=INK, ls="none", ms=8, mew=1.4),
          Line2D([0], [0], marker="o", mfc="white", mec=RED, mew=1.0, ls="none", ms=4.5),
          Line2D([0], [0], marker="o", mfc="white", mec=MUTED, mew=1.0, ls="none", ms=4.5)]
    ax.legend(hs, ["Pathogenic", "gnomAD", "gnomAD, seen twice or more",
                   "Pathogenic, outside low complexity", "gnomAD, outside low complexity"], ncol=2,
              loc="upper center", bbox_to_anchor=(0.30, -0.36), fontsize=6.5,
              handletextpad=0.2, columnspacing=0.8)
    tidy(ax, xlab="Frameshift share", title="c  Recurrence and low complexity", grid="x")
    save(fig, "figure4")


def ecdf(ax, vals, color, label=None, lw=1.6, ls="-"):
    v = np.sort(np.array(vals, dtype=float))
    y = np.arange(1, len(v) + 1) / len(v)
    ax.step(v, y, where="post", color=color, lw=lw, label=label, ls=ls)


def figure5():
    tails = defaultdict(lambda: {"p": [], "g": []})
    sy = {"p": [], "g": []}
    for r in VT:
        if r["variant_class"] != "frameshift" or r["fsd"] == "":
            continue
        k = "p" if r["source"] == "clinvar_pathogenic" else "g"
        tails[r["gene"]][k].append(int(r["fsd"]))
        if r["gene"] == "SYNGAP1" and r["codon"]:
            sy[k].append((int(r["codon"]), int(r["fsd"])))
    genes = W["frame"]["neo_tail_q_lt_0.05"]["mwu"]
    fq = load("derived_fsd_q.json")["mwu"]["genes"]
    fig = plt.figure(figsize=(7.4, 5.0))
    gsp = fig.add_gridspec(2, len(genes), height_ratios=[1.0, 1.3], hspace=0.85, wspace=0.6)
    for i, g in enumerate(genes):
        ax = fig.add_subplot(gsp[0, i])
        assert len(tails[g]["p"]) == F2[g]["n_pathogenic"], g
        assert len(tails[g]["g"]) == F2[g]["n_population"], g
        ecdf(ax, tails[g]["g"], MUTED)
        ecdf(ax, tails[g]["p"], RED)
        ax.set_xscale("symlog", linthresh=10)
        ax.set_ylim(0, 1.02)
        ax.set_xlim(0, max(max(tails[g]["p"]), max(tails[g]["g"])) * 1.15)
        ax.set_xticks([0, 10, 100])
        ax.set_xticklabels(["0", "10", "100"])
        tidy(ax, ylab="Cumulative share" if i == 0 else None,
             title=("a  " if i == 0 else "") + g, grid="y")

    ax = fig.add_subplot(gsp[1, :3])
    carry = W["pooled"]["neo_tail"]["carrying_gene"]
    allp, allg, xp, xg = [], [], [], []
    for g, d in tails.items():
        if g not in ISO:
            continue
        allp += d["p"]; allg += d["g"]
        if g != carry:
            xp += d["p"]; xg += d["g"]
    pa = W["pooled"]["neo_tail"]["pooled_all"]
    assert (len(allp), len(allg)) == (pa["n_path"], pa["n_pop"]), (len(allp), len(allg))
    ecdf(ax, allg, MUTED, "gnomAD")
    ecdf(ax, allp, RED, "Pathogenic")
    ecdf(ax, xg, MUTED, f"gnomAD without {carry}", lw=1.0, ls="--")
    ecdf(ax, xp, RED, f"Pathogenic without {carry}", lw=1.0, ls="--")
    ax.set_xscale("symlog", linthresh=10)
    top = max(max(allp), max(allg)) * 1.15
    ax.set_xlim(0, top)
    ticks = [t for t in (0, 10, 100, 1000) if t <= top]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower right", fontsize=6.5, handlelength=1.6)
    tidy(ax, xlab="Novel residues before the new stop", ylab="Cumulative share",
         title="b  Pooled over the isolated tier", grid="y")

    ax = fig.add_subplot(gsp[1, 3:])
    esc = FR["SYNGAP1"]["nmd_escape_from_codon"]
    L = FR["SYNGAP1"]["protein_aa"]
    ax.axvspan(esc, L, color="#e9e8e4", lw=0)
    for k, col in (("g", MUTED), ("p", RED)):
        ax.scatter([a for a, _ in sy[k]], [b for _, b in sy[k]], s=7, color=col,
                   edgecolor="none", alpha=0.85)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_yticks([0, 10, 100])
    ax.set_yticklabels(["0", "10", "100"])
    ax.set_xlim(0, L * 1.02)
    tidy(ax, xlab="SYNGAP1 codon", ylab="Novel residues", title="c  SYNGAP1", grid="y")
    save(fig, "figure5")


def strip_credentials():
    import struct
    for f in sorted(os.listdir(FIGS)):
        p = os.path.join(FIGS, f)
        if f.endswith(".png"):
            b = open(p, "rb").read()
            out, i = [b[:8]], 8
            while i < len(b):
                n = struct.unpack(">I", b[i:i + 4])[0]
                if b[i + 4:i + 8] not in (b"caBX", b"jumb"):
                    out.append(b[i:i + 12 + n])
                i += 12 + n
            open(p, "wb").write(b"".join(out))
        elif f.endswith(".svg"):
            s = open(p, encoding="utf-8").read()
            s = re.sub(r"<c2pa:manifest>.*?</c2pa:manifest>", "", s, flags=re.S)
            s = s.replace(' xmlns:c2pa="http://c2pa.org/manifest"', "")
            open(p, "w", encoding="utf-8").write(s)
        else:
            continue
        left = open(p, "rb").read().lower()
        if any(m in left for m in (b"c2pa", b"jumbf", b"cabx", b"contentcredentials")):
            raise SystemExit(f"content credentials remain in {f}")
    print("content credentials stripped from every figure")


def record_environment():
    import platform
    import matplotlib
    env = {"python": platform.python_version(), "numpy": np.__version__,
           "matplotlib": matplotlib.__version__}
    with open(os.path.join(DATA, "figure_environment.json"), "w") as fh:
        json.dump(env, fh, indent=1)


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    strip_credentials()
    record_environment()
