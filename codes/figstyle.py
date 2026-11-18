import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "..", "figures")

S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
RED = "#d62728"
GREY = "#b9b8b3"
WHITE = "#ffffff"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"

TITLE_SIZE = 9.5
TITLE_PAD = 7

plt.rcParams.update({
    "figure.facecolor": WHITE, "axes.facecolor": WHITE, "savefig.facecolor": WHITE,
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.7,
    "axes.labelcolor": INK, "axes.titlesize": TITLE_SIZE,
    "axes.titleweight": "bold",
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "grid.color": "#e9e8e4", "grid.linewidth": 0.5,
    "lines.linewidth": 2.0, "lines.markersize": 4.5,
})


def tidy(ax, ylab=None, xlab=None, title=None, grid="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if grid:
        ax.set_axisbelow(True)
        ax.grid(axis=grid, which="major")
    if ylab:
        ax.set_ylabel(ylab)
    if xlab:
        ax.set_xlabel(xlab)
    if title:
        ax.set_title(title, loc="left", color=INK, pad=TITLE_PAD)


def panel_extent(ax, renderer):

    from matplotlib.transforms import Bbox
    boxes = [ax.get_tightbbox(renderer)]
    for lbl in (ax.yaxis.label, ax.xaxis.label):
        if lbl.get_text().strip():
            boxes.append(lbl.get_window_extent(renderer))
    for axis in ("y", "x"):
        for lbl in drawn_ticklabels(ax, axis):
            boxes.append(lbl.get_window_extent(renderer))
    return Bbox.union([b for b in boxes if b.width and b.height])


def panel_headings(fig):

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    gap = (TITLE_PAD / 72.0) / fig.get_figheight()
    moved = []
    for ax in fig.axes:
        t = ax.get_title(loc="left")
        if not t:
            continue
        ax.set_title("", loc="left")
        moved.append((ax, t))
    if not moved:
        return


    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    box = {ax: inv.transform_bbox(panel_extent(ax, r)) for ax, _ in moved}


    rows_ = {}
    for ax, _ in moved:
        key = round(ax.get_position().y1, 3)
        rows_[key] = max(rows_.get(key, 0.0), box[ax].y1)


    cols = {}
    for ax, _ in moved:
        key = round(ax.get_position().x0, 3)
        cols[key] = min(cols.get(key, 1.0), box[ax].x0)


    left_of = {}
    for key in rows_:
        row = sorted((a for a, _ in moved
                      if round(a.get_position().y1, 3) == key),
                     key=lambda a: a.get_position().x0)
        for prev, ax in zip(row, row[1:]):
            left_of[ax] = box[prev].x1
    placed = {}
    for ax, t in moved:
        y = rows_[round(ax.get_position().y1, 3)] + gap
        x = cols[round(ax.get_position().x0, 3)]
        if ax in left_of:
            x = max(x, (left_of[ax] + ax.get_position().x0) / 2.0)
        placed[ax] = fig.text(x, y, t, ha="left", va="bottom",
                              color=INK, fontsize=TITLE_SIZE, fontweight="bold")


    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pad = (6.0 / 72.0) / fig.get_figwidth()
    for key in rows_:
        row = sorted((a for a, _ in moved
                      if round(a.get_position().y1, 3) == key),
                     key=lambda a: a.get_position().x0)
        right = None
        for ax in row:
            txt = placed[ax]
            ext = inv.transform_bbox(txt.get_window_extent(r))
            if right is not None and ext.x0 < right + pad:


                txt.set_x(ax.get_position().x0)
                fig.canvas.draw()
                ext = inv.transform_bbox(txt.get_window_extent(r))
                if ext.x0 < right + pad:
                    txt.set_x(right + pad)
                    fig.canvas.draw()
                    ext = inv.transform_bbox(txt.get_window_extent(r))
            right = ext.x1


def drawn_ticklabels(ax, axis):

    lo, hi = ax.get_ylim() if axis == "y" else ax.get_xlim()
    lo, hi = min(lo, hi), max(lo, hi)
    locs = ax.get_yticks() if axis == "y" else ax.get_xticks()
    labs = ax.get_yticklabels() if axis == "y" else ax.get_xticklabels()
    return [l for loc, l in zip(locs, labs)
            if lo - 1e-9 <= loc <= hi + 1e-9 and l.get_text().strip()]


def check_headings(fig, name):

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    heads = [(t.get_text(), t.get_window_extent(r)) for t in fig.texts]
    others = []
    for ax in fig.axes:
        items = (list(ax.texts) + drawn_ticklabels(ax, "y")
                 + drawn_ticklabels(ax, "x"))
        items += [l for l in (ax.yaxis.label, ax.xaxis.label)
                  if l.get_text().strip()]
        lg = ax.get_legend()
        if lg is not None:
            items += list(lg.get_texts())
        others += [(t.get_text(), t.get_window_extent(r))
                   for t in items if t.get_text().strip()]
    bad = []
    for hs, he in heads:
        bad += ["%s: heading %r lands on %r" % (name, hs, o)
                for o, oe in others if he.overlaps(oe)]
    for i in range(len(heads)):
        for j in range(i + 1, len(heads)):
            if heads[i][1].overlaps(heads[j][1]):
                bad.append("%s: heading %r lands on heading %r"
                           % (name, heads[i][0], heads[j][0]))


    for ax in fig.axes:
        lab = [t for t in ax.texts if t.get_text().strip()]
        for i in range(len(lab)):
            for j in range(i + 1, len(lab)):
                if lab[i].get_window_extent(r).overlaps(
                        lab[j].get_window_extent(r)):
                    bad.append("%s: label %r lands on label %r"
                               % (name, lab[i].get_text(), lab[j].get_text()))


    from matplotlib.transforms import Bbox
    drawn = {}
    for n, ax in enumerate(fig.axes):
        bs = []
        for art in list(ax.lines) + list(ax.collections) + list(ax.patches):
            if not art.get_visible():
                continue
            try:
                e = art.get_window_extent(r)
            except Exception:
                continue


            e = Bbox.intersection(e, ax.get_window_extent(r))
            if e is not None and e.width and e.height:
                bs.append(e)
        drawn[n] = bs
    for n, ax in enumerate(fig.axes):
        items = ([t for t in ax.texts if t.get_text().strip()]
                 + drawn_ticklabels(ax, "y") + drawn_ticklabels(ax, "x"))
        for t in items:
            e = t.get_window_extent(r)
            for m, marks in drawn.items():
                if m == n:
                    continue
                if any(e.overlaps(b) for b in marks):
                    bad.append("%s: %r from one panel lands on another panel's "
                               "drawing" % (name, t.get_text()))
                    break
    if bad:
        raise SystemExit("\n  ".join(["fatal: overlapping text in %s" % name] + bad))


def save(fig, name):
    panel_headings(fig)
    check_headings(fig, name)
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(os.path.join(FIGS, name + ".png"), dpi=400, bbox_inches="tight",
                facecolor=WHITE)
    fig.savefig(os.path.join(FIGS, name + ".svg"), bbox_inches="tight",
                facecolor=WHITE)
    plt.close(fig)
    print("wrote", name + ".png and .svg")
