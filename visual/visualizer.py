#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")   # non-GUI backend

import matplotlib.pyplot as plt
import re, sys, os
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ─── CONFIG: fixed grid size ─────────────────────────────────────────────────
M, N = 4, 5   # rows, columns

# ─── ARGS & PATHS ─────────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: python3 visualizer.py <filename_without_extension>")
    sys.exit(1)

base     = sys.argv[1]
in_path  = f"../output/{base}.out"
out_path = f"../visual/{base}.png"

if not os.path.exists(in_path):
    print(f"Input file not found: {in_path}")
    sys.exit(1)

# ─── READ & PARSE ONE MATRIX PAIR ─────────────────────────────────────────────
horiz = []
vert  = []
state = 0   # 0=seek H tag, 1=read H, 2=seek V tag, 3=read V

with open(in_path) as f:
    for raw in f:
        line = raw.strip()
        if not line:
            continue

        if state == 0 and line.lower().startswith("horizontal matrix"):
            state = 1
            continue

        if state == 1:
            if re.fullmatch(r"[0-9 ]+", line):
                horiz.append(list(map(int, line.split())))
                if len(horiz) == M:
                    state = 2
            continue

        if state == 2 and line.lower().startswith("vertical matrix"):
            state = 3
            continue

        if state == 3:
            if re.fullmatch(r"[0-9 ]+", line):
                vert.append(list(map(int, line.split())))
                if len(vert) == M - 1:
                    break
            continue

# sanity check
if len(horiz) != M or len(vert) != M-1:
    sys.exit(f"Parsing error: got {len(horiz)}×{len(horiz[0])} horiz and {len(vert)}×{len(vert[0])} vert, expected {M}×{N-1} and {M-1}×{N}")

# ─── PLOTTING ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(N, M))
ax.set_xlim(-0.5, N-0.5)
ax.set_ylim(-0.5, M-0.5)
ax.set_aspect("equal")

# draw grid points
for y in range(M):
    for x in range(N):
        ax.plot(x, y, 'ko', ms=4)

# fixed colormap ranges
cmaph, cmapv = plt.cm.Greens, plt.cm.Reds
norm_h = Normalize(vmin=0, vmax=4, clip=True)  # always 0–4
norm_v = Normalize(vmin=0, vmax=3, clip=True)  # always 0–3

# draw horizontal edges
for y in range(M):
    for x in range(N-1):
        c = horiz[y][x]
        ax.plot([x, x+1], [y, y],
                color=cmaph(norm_h(c)), lw=2)
        ax.text(x+0.5, y-0.2, str(c),
                ha='center', va='center', fontsize=8)

# draw vertical edges
for y in range(M-1):
    for x in range(N):
        c = vert[y][x]
        ax.plot([x, x], [y, y+1],
                color=cmapv(norm_v(c)), lw=2)
        ax.text(x+0.1, y+0.5, str(c),
                ha='left', va='center', fontsize=8)

ax.axis('off')

# two small colourbars, matching the fixed ranges
sm_h = plt.cm.ScalarMappable(cmap=cmaph, norm=norm_h)
sm_h.set_array([])
sm_v = plt.cm.ScalarMappable(cmap=cmapv, norm=norm_v)
sm_v.set_array([])

divider = make_axes_locatable(ax)
cax1 = divider.append_axes("right", size="5%", pad=0.05)
cax2 = divider.append_axes("right", size="5%", pad=0.7)
plt.colorbar(sm_h, cax=cax1, label="Horizontal (0–4)")
plt.colorbar(sm_v, cax=cax2, label="Vertical (0–3)")

plt.title(f"Routing Grid ({M}×{N})")
plt.tight_layout()

os.makedirs(os.path.dirname(out_path), exist_ok=True)
plt.savefig(out_path)
print(f"Saved image to {out_path}")
