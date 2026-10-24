#!/usr/bin/env python3
"""Draws the twenty points of the grid example (Supplementary Figure S6a): the points of
examples/manhattan_distances_20x20-points.tsv, whose Manhattan distances are examples/manhattan_distances_20x20.dist,
with their taxon labels. The points were matched to the taxa by the distances (the matching is unique).
Writes manuscript/figures/manhattan_distances_20x20-points.pdf. Needs matplotlib.
  python3 manuscript/experiments/manhattan-points.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
POINTS = os.path.join(ROOT, "examples", "manhattan_distances_20x20-points.tsv")
OUT = os.path.join(ROOT, "manuscript", "figures", "manhattan_distances_20x20-points.pdf")

# where the label of a point goes (dx, dy, horizontal alignment); the default is to the upper right
PLACEMENT = {"t4": (-0.45, 0, "right"), "t5": (0, -0.95, "center"), "t9": (-0.45, 0.2, "right")}

points = []
with open(POINTS) as f:
    next(f)
    for line in f:
        taxon, x, y = line.split()
        points.append((taxon, int(x), int(y)))

fig, ax = plt.subplots(figsize=(4.2, 4.2))
ax.set_xticks(range(0, 21), minor=True)
ax.set_yticks(range(0, 21), minor=True)
ax.grid(True, which="both", color="0.85", linewidth=0.5)
ax.set_axisbelow(True)
ax.scatter([x for _, x, _ in points], [y for _, _, y in points], s=28, color="black", zorder=3)
for taxon, x, y in points:
    dx, dy, ha = PLACEMENT.get(taxon, (0.4, 0.3, "left"))
    ax.text(x + dx, y + dy, taxon, fontsize=10, ha=ha, va="center", zorder=4)
ax.set_xlim(-1, 21.6)
ax.set_ylim(-1, 21.2)
ax.set_xticks(range(0, 21, 5))
ax.set_yticks(range(0, 21, 5))
ax.tick_params(labelsize=8)
ax.set_aspect("equal")
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
fig.tight_layout(pad=0.3)
fig.savefig(OUT)
print("written", OUT)
