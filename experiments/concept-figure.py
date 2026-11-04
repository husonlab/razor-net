#!/usr/bin/env python3
"""Draws panels (b) to (d) of Supplementary Figure S7, the idea behind RazorNet, on the five-point example
examples/five.dist (panel (a) is its distance matrix, set in LaTeX): (b) the taxa-only realization, (c) the
compactification of the slack point a, (d) the network of length 8 that RazorNet returns. Positions are placed
by hand; the numbers come from the metric. Writes experiments/figures/concept-{b,c,d}.pdf. Needs matplotlib.
  python3 experiments/concept-figure.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURES = os.path.join(HERE, "figures")
TAXON = dict(s=70, color="black", zorder=4)
AUX = dict(s=38, facecolors="white", edgecolors="black", linewidths=1.0, zorder=4)
LABEL = dict(fontsize=11, ha="center", va="center", zorder=5)
WEIGHT = dict(fontsize=9, ha="center", va="center", color="0.3", zorder=5,
              bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none"))


def panel(name, nodes, edges, labels, size=(3.6, 3.0), xlim=None, ylim=None):
    fig, ax = plt.subplots(figsize=size)
    for edge in edges:
        u, v, w, style = edge[:4]
        t = edge[4] if len(edge) > 4 else 0.5   # where along the edge the weight is written
        (x1, y1), (x2, y2) = nodes[u], nodes[v]
        if style == "solid":
            ax.plot([x1, x2], [y1, y2], color="black", linewidth=1.1, zorder=2)
        elif style == "faint":
            ax.plot([x1, x2], [y1, y2], color="0.7", linewidth=0.9, linestyle=(0, (3, 2)), zorder=1)
        if w is not None:
            ax.text(x1 + t * (x2 - x1), y1 + t * (y2 - y1), str(w), **WEIGHT)
    for n, (x, y) in nodes.items():
        if n.startswith("u"):
            ax.scatter([x], [y], **AUX)
        else:
            ax.scatter([x], [y], **TAXON)
    for n, (dx, dy, text) in labels.items():
        x, y = nodes[n]
        ax.text(x + dx, y + dy, text, **LABEL)
    ax.set_aspect("equal")
    ax.axis("off")
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    fig.tight_layout(pad=0.2)
    fig.savefig(os.path.join(FIGURES, "concept-%s.pdf" % name))
    plt.close(fig)


# (b) the taxa-only realization: every pair with no third taxon between it, weighted by its distance; a and d
# are not joined because e lies between them. Length 24.
nodes = {"a": (0, 0), "b": (0.9, 2.5), "c": (3.1, 2.5), "d": (4, 0), "e": (2, -1.6)}
edges = [("a", "b", 3, "solid"), ("a", "c", 3, "solid"), ("a", "e", 2, "solid"), ("b", "c", 2, "solid"),
         ("b", "d", 3, "solid"), ("b", "e", 3, "solid"), ("c", "d", 3, "solid"), ("c", "e", 3, "solid"),
         ("d", "e", 2, "solid")]
labels = {"a": (-0.4, 0, "a"), "b": (-0.35, 0.3, "b"), "c": (0.35, 0.3, "c"), "d": (0.4, 0, "d"), "e": (0, -0.45, "e")}
panel("b", nodes, edges, labels, xlim=(-0.8, 4.8), ylim=(-2.3, 3.1))

# (c) the compactification of a: its slack is 1, attained by b and e; the new point a' is the median of a, b
# and e, at distance 1 from a, 2 from b and 1 from e; its distances to c and d follow (2 and 3).
nodes = {"a": (0, 0), "ua": (1.0, 0), "b": (2.0, 1.7), "e": (1.7, -1.2), "c": (4.3, 1.2), "d": (4.3, -1.2)}
edges = [("a", "b", 3, "faint", 0.5), ("a", "e", 2, "faint", 0.5), ("b", "e", 3, "faint", 0.82),
         ("a", "ua", 1, "solid"), ("ua", "b", 2, "solid"), ("ua", "e", 1, "solid"),
         ("ua", "c", 2, "faint", 0.72), ("ua", "d", 3, "faint", 0.72)]
labels = {"a": (-0.4, 0, "a"), "ua": (0.05, 0.42, "a′"), "b": (0.0, 0.42, "b"), "e": (0.0, -0.45, "e"),
          "c": (0.4, 0, "c"), "d": (0.4, 0, "d")}
panel("c", nodes, edges, labels, xlim=(-0.8, 5.0), ylim=(-1.9, 2.3))

# (d) the network RazorNet returns: length 8, three additional nodes, every edge of length 1, one cycle
nodes = {"a": (0, 0), "u7": (1, 0), "e": (2, 0), "u8": (3, 0), "d": (4, 0), "u6": (2, 1.4), "b": (1, 2.4), "c": (3, 2.4)}
edges = [("a", "u7", 1, "solid"), ("u7", "e", 1, "solid"), ("e", "u8", 1, "solid"), ("u8", "d", 1, "solid"),
         ("u6", "u7", 1, "solid"), ("u6", "u8", 1, "solid"), ("b", "u6", 1, "solid"), ("c", "u6", 1, "solid")]
labels = {"a": (-0.4, 0, "a"), "b": (-0.35, 0.3, "b"), "c": (0.35, 0.3, "c"), "d": (0.4, 0, "d"), "e": (0, -0.45, "e")}
panel("d", nodes, edges, labels, xlim=(-0.8, 4.8), ylim=(-0.8, 3.0))
print("written concept-b/c/d.pdf")
