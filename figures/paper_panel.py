"""One panel in the style of the paper figures (axes x, y, patch numbers, no title).

usage: py paper_panel.py <out.png> <prefix.txt>
"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from draw_lines import read_lines

out, path = sys.argv[1], sys.argv[2]
lines = read_lines(path)
fig, ax = plt.subplots(figsize=(7, 7))
cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
pts = {}
for line in lines:
    p = line[0][2]
    ax.plot([q[0] for q in line], [q[1] for q in line], color=cols[p % len(cols)], lw=0.6)
    pts.setdefault(p, []).extend((q[0], q[1]) for q in line)
for p, P in pts.items():
    c = np.mean(np.array(P), axis=0)
    ax.text(c[0], c[1], str(p), ha='center', va='center', fontsize=11, fontweight='bold')
ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')
fig.tight_layout()
fig.savefig(out, dpi=200, bbox_inches='tight')
