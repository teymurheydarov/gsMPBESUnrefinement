r"""Figure 8 of the paper (fig:rotor_subdomains) from the published module: the first candidate
group and its fitting subdomain Delta(c) for lambda = 0, lambda = 2 and the global setting.

Sources (nothing is recomputed):
  data/first_candidate_<run>.txt        written by unrefinement_example --dump_first (the runs
                                        rotor_*_ef00002 with --dump_first first_candidate.txt):
                                        the cells of the first group and the boxes of Delta(c)
                                        (union of the supports of the free functions) per patch
  data/rotorRE_lines.txt, data/rotorRE_lines_grid.txt: patch outlines and 64x64 grids of
                                        data/geometries/rotorRE.xml, written by the research code
  rotorRE.xml: patch 22 is on level 4 (0-based) except the box (0.25, 0.75)^2 on level 5, all other
                                        patches on level 3

Run:  py figures\plot_subdomains.py
Output: figures/output/rotor_subdomains_{geometry,lambda0,lambda2,global}.png
"""
import io
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
LINES = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "output")
FREE, FROZEN, EDGE, CELL, ONE, TWO, COARSE = "#3b6ea5", "#ececec", "#7a7a7a", "#b3271e", "#f2c9a0", "#d9822b", "#dce9f5"
FINE_PATCH, FINE_BOX = 22, (0.25, 0.75, 0.25, 0.75)   # (u0, u1, v0, v1) of the level-5 box of patch 22
PANELS = [("lambda0_ef0002", "lambda0", r"$\lambda=0$"),
          ("lambda2_ef0002", "lambda2", r"$\lambda=2$"),
          ("global_ef0002", "global", "global")]


def outlines():
    out, cur = [], None
    for line in io.open(os.path.join(LINES, "rotorRE_lines.txt"), encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("patch"):
            cur = []
            out.append(cur)
            continue
        x, y = line.split()
        cur.append((float(x), float(y)))
    return out


def grids():
    g, cur = [], None
    for line in io.open(os.path.join(LINES, "rotorRE_lines_grid.txt"), encoding="utf-8"):
        line = line.strip()
        if line.startswith("patch"):
            cur = {"n": int(line.split()[2]), "pts": []}
            g.append(cur)
            continue
        if line:
            x, y = line.split()
            cur["pts"].append((float(x), float(y)))
    return g


def to_phys(grid, u, v):
    """Bilinear interpolation of the patch map on its (n+1)^2 grid."""
    n = grid["n"]
    s, t = min(max(u, 0.0), 1.0) * n, min(max(v, 0.0), 1.0) * n
    i, j = min(int(s), n - 1), min(int(t), n - 1)
    a, b = s - i, t - j
    P = grid["pts"]
    p00, p10 = P[j * (n + 1) + i], P[j * (n + 1) + i + 1]
    p01, p11 = P[(j + 1) * (n + 1) + i], P[(j + 1) * (n + 1) + i + 1]
    return ((1 - a) * (1 - b) * p00[0] + a * (1 - b) * p10[0] + (1 - a) * b * p01[0] + a * b * p11[0],
            (1 - a) * (1 - b) * p00[1] + a * (1 - b) * p10[1] + (1 - a) * b * p01[1] + a * b * p11[1])


def box_polygon(grid, box, k=16):
    u0, u1, v0, v1 = box
    pts = []
    for e in range(4):
        for m in range(k):
            t = m / float(k)
            u, v = [(u0 + t * (u1 - u0), v0), (u1, v0 + t * (v1 - v0)),
                    (u1 - t * (u1 - u0), v1), (u0, v1 - t * (v1 - v0))][e]
            pts.append(to_phys(grid, u, v))
    return pts


def first_candidate(run):
    d = {"cells": [], "domain": [], "free": None}
    for line in io.open(os.path.join(LINES, "first_candidate_" + run + ".txt"), encoding="utf-8"):
        w = line.split()
        if not w or w[0].startswith("#"):
            continue
        if w[0] == "free":
            d["free"] = (int(w[1]), int(w[3]))
        elif w[0] in ("cell", "domain"):
            d["cells" if w[0] == "cell" else "domain"].append((int(w[1]), tuple(float(x) for x in w[2:6])))
    return d


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "rotor_subdomains_" + name + ".png")
    fig.savefig(path, dpi=165, bbox_inches="tight")
    plt.close(fig)
    print("written:", path)


ps, gs = outlines(), grids()
data = {run: first_candidate(run) for run, _, _ in PANELS}
for run, _, _ in PANELS:
    d = data[run]
    print("%-16s free %s, %d cells on patches %s, %d domain boxes on patches %s" % (
        run, d["free"], len(d["cells"]), sorted({p for p, _ in d["cells"]}),
        len(d["domain"]), sorted({p for p, _ in d["domain"]})))
groups = {tuple(sorted(data[run]["cells"])) for run, _, _ in PANELS}
assert len(groups) == 1, "the first group differs between the settings"
cells = data[PANELS[0][0]]["cells"]

# zoom window: the patches of the largest local subdomain (lambda = 2), padded
zoom = sorted({p for p, _ in data["lambda2_ef0002"]["domain"]}) or [FINE_PATCH]
xs = [x for p in zoom for x, _ in ps[p]]
ys = [y for p in zoom for _, y in ps[p]]
cx, cy = 0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys))
w = 1.35 * max(max(xs) - min(xs), max(ys) - min(ys))
win = (cx - w, cx + w, cy - w, cy + w)

# geometry: patch 22 light (one additional level), its level-5 box dark (two), zoom window dashed
fig, ax = plt.subplots(figsize=(4.5, 4.9))
for i, p in enumerate(ps):
    ax.add_patch(Polygon(p, closed=True, facecolor=ONE if i == FINE_PATCH else COARSE, edgecolor=EDGE, linewidth=0.4))
ax.add_patch(Polygon(box_polygon(gs[FINE_PATCH], FINE_BOX), closed=True, facecolor=TWO, edgecolor="none"))
ax.add_patch(Rectangle((win[0], win[2]), win[1] - win[0], win[3] - win[2], fill=False, edgecolor="k",
                       linewidth=0.8, linestyle="--"))
ax.set_title("rotor", fontsize=12)
ax.relim()
ax.autoscale_view()
ax.set_aspect("equal")
ax.axis("off")
save(fig, "geometry")

for run, name, title in PANELS:
    d = data[run]
    local = bool(d["domain"])
    fig, ax = plt.subplots(figsize=(4.5, 4.9))
    for p in ps:
        ax.add_patch(Polygon(p, closed=True, facecolor=FROZEN if local else FREE, edgecolor="none"))
    for patch, box in d["domain"]:
        ax.add_patch(Polygon(box_polygon(gs[patch], box), closed=True, facecolor=FREE, edgecolor="none"))
    for p in ps:
        ax.add_patch(Polygon(p, closed=True, facecolor="none", edgecolor=EDGE, linewidth=0.4))
    for patch, box in cells:
        ax.add_patch(Polygon(box_polygon(gs[patch], box), closed=True, facecolor="none", edgecolor=CELL, linewidth=0.9))
    ax.set_title(title, fontsize=12)
    ax.set_xlim(win[0], win[1])
    ax.set_ylim(win[2], win[3])
    ax.set_aspect("equal")
    ax.axis("off")
    save(fig, name)
