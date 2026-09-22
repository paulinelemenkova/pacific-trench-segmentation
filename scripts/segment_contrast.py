#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator

CSV = sys.argv[1] if len(sys.argv) > 1 else "boundary_contrast.csv"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Arial", "Helvetica",
                        "Liberation Sans", "TeX Gyre Heros"],
    "font.size": 9,
    "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 9,
    "axes.linewidth": 0.7, "xtick.direction": "out", "ytick.direction": "out",
    "figure.dpi": 150, "savefig.dpi": 600,
})

SECTOR = {
    "NW": ("Northwest Pacific", "#0072B2", "o"),
    "W":  ("West Pacific",      "#E69F00", "s"),
    "SW": ("Southwest Pacific", "#009E73", "^"),
    "E":  ("East Pacific",      "#D55E00", "D"),
}
GRID = "#B9BDC4"

VARS = [
    ("depth", "(a)", "Axial depth (km)",        r"$\Delta$ depth (km)"),
    ("width", "(b)", "Trench width (km)",       r"$\Delta$ width (km)"),
    ("asym",  "(c)", "Wall asymmetry",          r"$\Delta$ asymmetry"),
]

df = pd.read_csv(CSV, comment="#")
order = ["NW", "W", "SW", "E"]

fig = plt.figure(figsize=(7.3, 3.5))
gs = GridSpec(2, 3, figure=fig, height_ratios=[3.0, 1.15],
              hspace=0.52, wspace=0.34,
              left=0.075, right=0.99, top=0.86, bottom=0.135)

for k, (stem, tag, label, dlabel) in enumerate(VARS):
    pre = df[f"{stem}_pre"].to_numpy() if f"{stem}_pre" in df else df[f"{stem}_pre_km"].to_numpy()
    post = df[f"{stem}_post"].to_numpy() if f"{stem}_post" in df else df[f"{stem}_post_km"].to_numpy()

    ax = fig.add_subplot(gs[0, k])
    lo = min(pre.min(), post.min()); hi = max(pre.max(), post.max())
    pad = 0.05 * (hi - lo)
    lim = (lo - pad, hi + pad)

    ax.plot(lim, lim, ls="--", lw=0.8, color="0.35", zorder=1)

    for s in order:
        m = df["sector"].to_numpy() == s
        if not m.any():
            continue
        _, col, mk = SECTOR[s]
        ax.scatter(pre[m], post[m], s=15, c=col, marker=mk,
                   edgecolors="white", linewidths=0.3, alpha=0.9, zorder=3)

    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal", "box")
    ax.set_xlabel(f"preceding segment\n{label}")
    if k == 0:
        ax.set_ylabel(f"following segment\n{label}")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(which="both", top=True, right=True)
    ax.annotate(tag, xy=(0.04, 0.93), xycoords="axes fraction",
                fontsize=11, fontweight="bold", va="top", ha="left")

    axh = fig.add_subplot(gs[1, k])
    d = post - pre
    r = np.nanmax(np.abs(d)); r = r if r > 0 else 1.0
    axh.hist(d, bins=22, range=(-r, r), color="0.55", edgecolor="white", lw=0.3)
    axh.axvline(0.0, color="0.15", lw=1.2, zorder=4)
    axh.set_xlim(-r, r)
    axh.set_xlabel(dlabel)
    if k == 0:
        axh.set_ylabel("count")
    axh.xaxis.set_minor_locator(AutoMinorLocator())
    axh.yaxis.set_minor_locator(AutoMinorLocator())
    axh.grid(True, which="major", axis="y", color=GRID, lw=0.5)
    axh.set_axisbelow(True)
    axh.tick_params(which="both", top=True, right=True)

handles = [Line2D([0], [0], ls="none", marker=SECTOR[s][2], markersize=6,
                  markerfacecolor=SECTOR[s][1], markeredgecolor="white",
                  markeredgewidth=0.3, label=SECTOR[s][0]) for s in order]
fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
           bbox_to_anchor=(0.5, 1.005), handletextpad=0.4, columnspacing=1.4)

fig.savefig("fig_segment_contrast.pdf", bbox_inches="tight")
fig.savefig("fig_segment_contrast.png", bbox_inches="tight")
print("wrote fig_segment_contrast.pdf and fig_segment_contrast.png; "
      f"n = {len(df)} boundaries")
