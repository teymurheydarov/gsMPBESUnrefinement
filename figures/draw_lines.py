"""Draw the line file of export_visualization_mesh (<prefix>.txt: "x y patch dir" per point,
blank line between lines) without the heuristics of viz_main.py.

usage: py draw_lines.py <out.png> "<title>" <prefix.txt> [<prefix.txt> ...]  (one panel per file;
       a panel title may be given as path=title)
"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_lines(path):
    lines, cur = [], []
    for raw in open(path, encoding='utf-8'):
        s = raw.strip()
        if not s or s.startswith('#'):
            if cur:
                lines.append(cur)
                cur = []
            continue
        x, y, p, d = s.split()
        cur.append((float(x), float(y), int(p), d))
    if cur:
        lines.append(cur)
    return lines


def main():
    out, title, specs = sys.argv[1], sys.argv[2], sys.argv[3:]
    fig, axs = plt.subplots(1, len(specs), figsize=(7 * len(specs), 7.5), squeeze=False)
    cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for ax, spec in zip(axs[0], specs):
        path, _, sub = spec.partition('=')
        for line in read_lines(path):
            p = line[0][2]
            ax.plot([q[0] for q in line], [q[1] for q in line], color=cols[p % len(cols)], lw=0.5)
        ax.set_aspect('equal')
        ax.set_title(sub or path.replace('\\', '/').split('/')[-1], fontsize=11)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    fig.savefig(out, dpi=110)


if __name__ == "__main__":
    main()
