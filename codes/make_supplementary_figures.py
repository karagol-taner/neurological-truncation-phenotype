#!/usr/bin/env python3

import json
import math
import os
import platform

import numpy as np
import matplotlib
from matplotlib.lines import Line2D

from figstyle import plt, tidy, save, RED, GREY, INK, INK2, MUTED

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
SUBUNITS = [("brain_epilepsy", "Epilepsy"), ("brain_neurodevelopmental", "Neurodevelopmental"),
            ("brain_neurodegeneration", "Neurodegeneration"), ("brain_malformation", "Malformation")]
RUNS = [("brain_whole", "All variants", INK, "o"),
        ("mane_whole", "MANE Select only", RED, "s"),
        ("lcrsym_whole", "Outside low complexity", INK2, "D"),
        ("manelcrsym_whole", "Both restrictions", MUTED, "^"),
        ("lcrflag_whole", "gnomAD flag, population only", None, "v")]


def load(*parts):
    with open(os.path.join(DATA, *parts)) as fh:
        return json.load(fh)


def run_file(run, name):
    return load(name) if run == "brain_whole" else load("sensitivity", run, name)


env = load("figure_environment.json")
here = {"python": platform.python_version(), "numpy": np.__version__, "matplotlib": matplotlib.__version__}
if env != here:
    raise SystemExit(f"figure environment {here} differs from the one the main figures used {env}")

RS = load("results_summary.json")


def seps_of(S):
    s_ = S["frame"]["q_lt_0.05"]
    return sorted(set(s_["T1"]) | set(s_["T4"]))


def tails_of(S):
    return S["frame"]["neo_tail_q_lt_0.05"]["mwu"]


def minq(e):
    qs = [e.get(k) for k in ("T1_q", "T4_q") if e and e.get(k) is not None]
    return min(qs) if qs else None


def dotplot(ax, genes, qof, title):
    thr = -math.log10(0.05)
    top = 0.0
    for i, g in enumerate(genes):
        y = len(genes) - 1 - i
        for j, (run, label, col, mk) in enumerate(RUNS):
            q = qof(run, g)
            if q is None:
                continue
            x = -math.log10(q)
            top = max(top, x)
            face, edge = ("white", INK) if col is None else (col, "none")
            ax.scatter([x], [y + (j - 2) * 0.13], s=16, marker=mk, facecolor=face, edgecolor=edge,
                       linewidths=0.8, zorder=3)
    ax.axvline(thr, color=MUTED, lw=0.7, ls="--")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(list(reversed(genes)))
    ax.set_ylim(-0.6, len(genes) - 0.4)
    ax.set_xlim(-0.05, top * 1.08)
    tidy(ax, xlab="−log₁₀ q", title=title, grid="x")


def figure_s1():
    frame = {run: run_file(run, "derived_frame.json")["genes"] for run, *_ in RUNS}
    tail = {run: run_file(run, "derived_fsd_q.json")["mwu"]["genes"] for run, *_ in RUNS}
    fgenes = sorted({g for run, *_ in RUNS for g in seps_of(RS[run])},
                    key=lambda g: (minq(frame["brain_whole"].get(g)) or 1.0, g))
    tgenes = sorted({g for run, *_ in RUNS for g in tails_of(RS[run])},
                    key=lambda g: (tail["brain_whole"].get(g, {}).get("q", 1.0), g))
    for run, *_ in RUNS:
        for g in seps_of(RS[run]):
            if minq(frame[run].get(g)) is None or minq(frame[run][g]) >= 0.05:
                raise SystemExit(f"{run}: {g} is listed as separating but its q does not pass")
        for g in tails_of(RS[run]):
            if tail[run][g]["q"] >= 0.05:
                raise SystemExit(f"{run}: {g} is listed as a tail gene but its q does not pass")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 4.4), gridspec_kw={"width_ratios": [1.15, 1]})
    dotplot(a, fgenes, lambda run, g: minq(frame[run].get(g)), "a  Reading frame")
    dotplot(b, tgenes, lambda run, g: tail[run][g]["q"] if g in tail[run] else None, "b  Novel-tail length")
    handles = [Line2D([], [], ls="", marker=mk, label=label, mfc="white" if col is None else col,
                      mec=INK if col is None else col, mew=0.8) for _, label, col, mk in RUNS]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.06))
    fig.subplots_adjust(bottom=0.2, wspace=0.45)
    save(fig, "figureS1")


def label_column(ax, items):

    if not items:
        return
    items = sorted(items, key=lambda t: -t[2])
    if len(items) == 1:
        g, x, y = items[0]
        ax.annotate(g, (x, y), xytext=(5, 4), textcoords="offset points", fontsize=6.5, color=INK)
        return
    step = 0.06
    xs = min(x for _, x, _ in items) - 0.16
    mid = sum(y for _, _, y in items) / len(items)
    y0 = min(0.97, mid + step * (len(items) - 1) / 2)
    for i, (g, x, y) in enumerate(items):
        yt = y0 - i * step
        ax.text(xs, yt, g, fontsize=6.5, color=INK, ha="right", va="center")
        ax.plot([xs + 0.01, x - 0.008], [yt, y], color=MUTED, lw=0.4, zorder=1)


def figure_s2():
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.8))
    jobs = []
    for ax, (u, label), letter in zip(axes.flat, SUBUNITS, "abcd"):
        G = load("units", u, "derived_frame.json")["genes"]
        genes = sorted(g for g, e in G.items() if e.get("T1_p") is not None)
        sep = seps_of(RS[u])
        if not set(sep) <= set(genes):
            raise SystemExit(f"{u}: a separating gene was not tested")
        for g in genes:
            e = G[g]
            on = g in sep
            col = (RED if e["path_fs_share"] < e["gnomad_fs_share"] else INK) if on else GREY
            ax.scatter([e["gnomad_fs_share"]], [e["path_fs_share"]], s=14 if on else 7, color=col,
                       edgecolor="none", alpha=0.95 if on else 0.5, zorder=3 if on else 2)
        ax.plot([0, 1], [0, 1], color=MUTED, lw=0.7, ls="--")
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.04)
        tidy(ax, xlab="gnomAD frameshift share", ylab="Pathogenic frameshift share",
             title=f"{letter}  {label}, {len(genes)} genes", grid="both")
        jobs.append((ax, [(g, G[g]["gnomad_fs_share"], G[g]["path_fs_share"]) for g in sep]))
    for ax, items in jobs:
        label_column(ax, items)
    fig.tight_layout()
    save(fig, "figureS2")


if __name__ == "__main__":
    figure_s1()
    figure_s2()
    from make_figures import strip_credentials
    strip_credentials()
