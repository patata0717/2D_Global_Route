import matplotlib
matplotlib.use("Agg")  # Use non-GUI backend

import matplotlib.pyplot as plt
import numpy as np
import re
import sys
import os
from matplotlib.colors import Normalize

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

# Clean empty lines
lines = [line.strip() for line in lines if line.strip()]

# Parse M and N
M = int(re.search(r"\d+", lines[0]).group())  # columns
N = int(re.search(r"\d+", lines[1]).group())  # rows

# Parse horizontal and vertical matrices
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
fig, ax = plt.subplots(figsize=(M, N))
ax.set_xlim(-0.5, M - 0.5)
ax.set_ylim(N - 0.5, - 0.5)
ax.set_aspect('equal')

# Flip Y-axis to make (0,0) bottom-left
flip_y = lambda y: (N - 1 - y)

# Normalize color by max edge count
max_val = max(max(map(max, horizontal)), max(map(max, vertical)))
cmap = plt.cm.hot

# Draw grid points
for y in range(N):
    for x in range(M):
        ax.plot(x, flip_y(y), 'ko', markersize=4)

# Draw horizontal edges
for y in range(N):
    for x in range(M - 1):
        count = horizontal[y][x]
        color = cmap(count / max_val)
        y_plot = flip_y(y)
        ax.plot([x, x + 1], [y_plot, y_plot], color=color, linewidth=2)
        ax.text(x + 0.5, y_plot - 0.1, str(count), ha='center', va='center', fontsize=8)

# Draw vertical edges
for y in range(N - 1):
    for x in range(M):
        count = vertical[y][x]
        color = cmap(count / max_val)
        y0 = flip_y(y)
        y1 = flip_y(y + 1)
        ax.plot([x, x], [y0, y1], color=color, linewidth=2)
        ax.text(x + 0.1, (y0 + y1) / 2, str(count), ha='left', va='center', fontsize=8)

# Hide axes
ax.axis('off')

# Add colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=0, vmax=max_val))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.7)
cbar.set_label("Edge Usage Count")

# Save the figure
plt.title(f"Superimposed Routing Grid ({M}×{N})")
plt.tight_layout()
plt.savefig(output_path)
print(f"Grid image saved to {output_path}")
