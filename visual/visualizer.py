import matplotlib
matplotlib.use("Agg")  # Use non-GUI backend

import matplotlib.pyplot as plt
import numpy as np
import re
import sys
import os
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

# --- Step 1: Get filename from argument ---
if len(sys.argv) != 2:
    print("Usage: python3 visualizer.py <filename_without_extension>")
    sys.exit(1)

filename = sys.argv[1]
input_path = f"../output/{filename}.out"
output_path = f"../visual/{filename}.png"

# --- Step 2: Read and parse the file ---
if not os.path.exists(input_path):
    print(f"Input file not found: {input_path}")
    sys.exit(1)

with open(input_path, "r") as f:
    lines = f.readlines()

lines = [line.strip() for line in lines if line.strip()]

M = int(re.search(r"\d+", lines[0]).group())  # rows
N = int(re.search(r"\d+", lines[1]).group())  # columns

horizontal = []
vertical = []
mode = None
for line in lines[2:]:
    if "Horizontal" in line:
        mode = "horizontal"
        continue
    elif "Vertical" in line:
        mode = "vertical"
        continue

    nums = list(map(int, line.split()))
    if mode == "horizontal":
        horizontal.append(nums)
    elif mode == "vertical":
        vertical.append(nums)

# --- Step 3: Plot the grid ---
fig, ax = plt.subplots(figsize=(N, M))
ax.set_xlim(-0.5, N - 0.5)
ax.set_ylim(-0.5, M - 0.5)  # (0,0) bottom-left
ax.set_aspect('equal')

# Draw grid points
for y in range(M):
    for x in range(N):
        ax.plot(x, y, 'ko', markersize=4)

# Choose a single colormap, but two different normalizations
cmaph = plt.cm.Greens
cmapv = plt.cm.Reds
norm_h = Normalize(vmin=0, vmax=4, clip=True)  # horizontal: 0–4
norm_v = Normalize(vmin=0, vmax=3, clip=True)  # vertical:   0–3

# Draw horizontal edges
for y in range(M):
    for x in range(N - 1):
        count = horizontal[y][x]
        color = cmaph(norm_h(count))
        ax.plot([x, x + 1], [y, y], color=color, linewidth=2)
        ax.text(x + 0.5, y - 0.2, str(count),
                ha='center', va='center', fontsize=8)

# Draw vertical edges
for y in range(M - 1):
    for x in range(N):
        count = vertical[y][x]
        color = cmapv(norm_v(count))
        ax.plot([x, x], [y, y + 1], color=color, linewidth=2)
        ax.text(x + 0.1, y + 0.5, str(count),
                ha='left', va='center', fontsize=8)

# Hide axes
ax.axis('off')

# --- Step 4: Add two colorbars ---
# Horizontal colorbar
sm_h = plt.cm.ScalarMappable(cmap=cmaph, norm=norm_h)
sm_h.set_array([])
# Vertical colorbar
sm_v = plt.cm.ScalarMappable(cmap=cmapv, norm=norm_v)
sm_v.set_array([])

divider = make_axes_locatable(ax)
cax_h = divider.append_axes("right", size="5%", pad=0.05)
cbar_h = plt.colorbar(sm_h, cax=cax_h)
cbar_h.set_label("Horizontal Usage Count (0–4)")

cax_v = divider.append_axes("right", size="5%", pad=0.7)
cbar_v = plt.colorbar(sm_v, cax=cax_v)
cbar_v.set_label("Vertical Usage Count (0–3)")

# Title & save
plt.title(f"Superimposed Routing Grid ({M}×{N})")
plt.tight_layout()
plt.savefig(output_path)
print(f"Grid image saved to {output_path}")
