#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import os

import matplotlib
import matplotlib.font_manager
import matplotlib.patches
import matplotlib.patheffects

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator

import analysis as an

HERE = os.path.dirname(os.path.abspath(__file__))
for _f in __import__("glob").glob(os.path.expanduser("~/.fonts/NimbusSans-*.otf")):
    matplotlib.font_manager.fontManager.addfont(_f)
OUT = os.path.join(HERE, "figures")
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Nimbus Sans", "Helvetica"],
        "font.size": 9,
        "axes.labelsize": 9.5,
        "axes.titlesize": 10,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.minor.size": 2.2,
        "ytick.minor.size": 2.2,
        "axes.labelpad": 2,
        "axes.titlepad": 3,
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.edgecolor": "0.6",
        "pdf.fonttype": 42,
        "savefig.dpi": 400,
    }
)
assert "dejavu" not in findfont(FontProperties(family=["Nimbus Sans"])).lower()

OI = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "verm": "#D55E00",
    "purple": "#CC79A7",
}
CLS = {
    "fracture_zone": ("Fracture zone", OI["blue"], "s"),
    "plateau": ("LIP plateau / ridge", OI["verm"], "D"),
    "seamount_ridge": ("Seamount", OI["green"], "^"),
}
SECTOR_COL = {
    "Northwest": "#FFFF00",
    "West": "#FF3030",
    "Southwest": "#7FFF00",
    "East": "#FF00FF",
}
SECTOR = {
    **{t: ("Northwest", SECTOR_COL["Northwest"]) for t in (1, 2, 3, 4, 5)},
    **{t: ("West", SECTOR_COL["West"]) for t in (6, 7, 8, 9, 10)},
    **{t: ("Southwest", SECTOR_COL["Southwest"]) for t in (11, 12, 13, 14, 15, 16, 17)},
    **{t: ("East", SECTOR_COL["East"]) for t in (18, 19, 20)},
}
EDGE = [
    matplotlib.patheffects.Stroke(linewidth=3.4, foreground="black"),
    matplotlib.patheffects.Normal(),
]
HITLAB = {
    "fracture_zone": "fracture zone",
    "plateau": "LIP plateau/ridge",
    "seamount_ridge": "seamount",
}
CREDIT = "Software: Python 3.11, NumPy, Matplotlib, Basemap. Source: authors."


def style(ax, grid="both"):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", top=True, right=True)
    if grid:
        ax.grid(
            True,
            which="major",
            axis=grid if grid != "both" else "both",
            lw=0.5,
            color="0.87",
        )
    ax.set_axisbelow(True)


def tag(ax, s, loc="left"):
    x = 0.015 if loc == "left" else 0.985
    ax.text(
        x,
        0.975,
        f"({s})",
        transform=ax.transAxes,
        ha=loc,
        va="top",
        fontsize=10.5,
        fontweight="bold",
        bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.85),
        zorder=20,
    )


def save(fig, name):
    fig.savefig(
        os.path.join(OUT, f"fig_{name}.pdf"), bbox_inches="tight", pad_inches=0.03
    )
    fig.savefig(
        os.path.join(OUT, f"fig_{name}.png"), bbox_inches="tight", pad_inches=0.03
    )
    plt.close(fig)


def rcsv(name, sub="results"):
    with open(os.path.join(HERE, sub, name)) as f:
        return list(csv.DictReader(f))


SUMMARY = {int(r["trench_id"]): r for r in rcsv("trench_summary.csv")}
CAT = rcsv("boundary_catalogue.csv")
CONTACTS = rcsv("contacts_in_windows.csv")
ALLC = rcsv("fabric_contacts.csv", "data")
COIN = rcsv("coincidence_by_class.csv")
SERIES = {t: an.load_series(t) for t in range(1, 21)}


def window(t):
    r = SUMMARY[t]
    if r["status"] != "analysed":
        return None
    s = SERIES[t]["s_km"]
    return int(np.searchsorted(s, float(r["window_start_km"]))), int(
        np.searchsorted(s, float(r["window_end_km"]))
    )


def basemap(
    ax,
    lon0=115,
    lon1=295,
    lat0=-56,
    lat1=63,
    relief="colour",
    res="l",
    dpar=20,
    dmer=30,
):
    from mpl_toolkits.basemap import Basemap

    m = Basemap(
        projection="mill",
        llcrnrlon=lon0,
        urcrnrlon=lon1,
        llcrnrlat=lat0,
        urcrnrlat=lat1,
        resolution=res,
        ax=ax,
    )
    if relief == "grey":
        import mpl_toolkits.basemap_data as bd
        from PIL import Image

        src = os.path.join(list(bd.__path__)[0], "etopo1.jpg")
        grey = os.path.join(OUT, "_etopo_grey.png")
        if not os.path.exists(grey):
            Image.open(src).convert("L").save(grey)
        m.warpimage(grey, scale=0.5, alpha=0.55, cmap="gray")
    else:
        m.etopo(scale=0.5, alpha=0.75)
    m.drawcoastlines(linewidth=0.3, color="0.25")
    kw = dict(linewidth=0.4, color="white", dashes=[1, 0], fontsize=8, zorder=1)
    m.drawparallels(np.arange(-60, 61, dpar), labels=[1, 0, 0, 0], **kw)
    m.drawmeridians(
        np.arange(120, 300, dmer),
        labels=[0, 0, 0, 1],
        fmt=lambda x: f"{x:.0f}°E" if x <= 180 else f"{360 - x:.0f}°W",
        **kw,
    )
    return m


def fig_studyarea():
    fig, ax = plt.subplots(figsize=(7.1, 4.6))
    m = basemap(ax)
    for t in range(1, 21):
        S = SERIES[t]
        w = window(t)
        lab, col = SECTOR[t]
        x, y = m(S["axis_lon"], S["axis_lat"])
        if w is None:
            m.plot(x, y, color="0.25", lw=1.2, ls=(0, (2, 1.5)), zorder=5)
        else:
            m.plot(x, y, color="0.35", lw=1.2, alpha=0.6, zorder=4)
            a, b = w
            m.plot(
                x[a : b + 1], y[a : b + 1], color=col, lw=2.0, zorder=5, path_effects=EDGE
            )
        i = len(x) // 2
        ax.annotate(
            str(t),
            (x[i], y[i]),
            xytext=(4, 3),
            textcoords="offset points",
            fontsize=8,
            fontweight="bold",
            color="black",
            zorder=6,
            path_effects=[
                matplotlib.patheffects.withStroke(linewidth=2, foreground="white")
            ],
        )
    names = [f"{t} {SERIES[t]['name']}" for t in range(1, 21)]
    h = [
        Line2D([], [], color=c, lw=2.0, label=s, path_effects=EDGE)
        for s, c in SECTOR_COL.items()
    ]
    h.append(
        Line2D(
            [], [], color="0.25", lw=1.2, ls=(0, (2, 1.5)), label="Excluded (trough)"
        )
    )
    leg = ax.legend(
        handles=h,
        loc="lower left",
        bbox_to_anchor=(0.46, 0.02),
        fontsize=8,
        title="Sector",
        title_fontsize=8.5,
    )
    leg.set_zorder(30)
    ax.text(
        1.01,
        1.0,
        "\n".join(names),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        linespacing=1.25,
    )
    save(fig, "studyarea")


def read_geometry():
    feats, cur, typ = [], [], None
    for ln in open(os.path.join(HERE, "data", "fabric_geometry.gmt")):
        ln = ln.strip()
        if ln.startswith(">"):
            if cur:
                feats.append((typ, np.array(cur)))
            typ = ln.split()[1]
            cur = []
        elif ln:
            cur.append(tuple(map(float, ln.split())))
    if cur:
        feats.append((typ, np.array(cur)))
    return feats


def fig_fabric():
    fig, ax = plt.subplots(figsize=(7.1, 4.6))
    m = basemap(ax, relief="grey")
    for typ, xy in read_geometry():
        x, y = m(xy[:, 0], xy[:, 1])
        if typ == "P":
            ax.fill(x, y, fc=OI["verm"], ec=OI["verm"], alpha=0.35, lw=0.8, zorder=3)
        else:
            m.plot(x, y, color=OI["blue"], lw=1.2, zorder=3)
    for t in range(1, 21):
        S = SERIES[t]
        x, y = m(S["axis_lon"], S["axis_lat"])
        m.plot(x, y, color="black", lw=1.2, zorder=4)
    for c in ("fracture_zone", "seamount_ridge", "plateau"):
        lab, col, mk = CLS[c]
        pts = [r for r in ALLC if r["class"] == c]
        x, y = m([float(r["lon"]) for r in pts], [float(r["lat"]) for r in pts])
        ax.scatter(
            x,
            y,
            s=22,
            marker=mk,
            c=col,
            edgecolors="black",
            linewidths=0.5,
            zorder=6,
            label=f"{lab} contact ({len(pts)})",
        )
    labs = {
        "Hikurangi Plateau": ((181.0, -41.5), (160.0, -33.0)),
        "Nazca Ridge": ((281.5, -17.5), (258.0, -11.0)),
        "Carnegie Ridge": ((276.5, -1.0), (258.0, 3.5)),
    }
    for nm, (xy, xyt) in labs.items():
        ax.annotate(
            nm,
            m(*xy),
            xytext=m(*xyt),
            fontsize=8.5,
            ha="center",
            va="center",
            zorder=7,
            arrowprops=dict(arrowstyle="-|>", lw=1.0, color="black"),
            path_effects=[
                matplotlib.patheffects.withStroke(linewidth=2, foreground="white")
            ],
        )
    h = [
        Line2D([], [], color="black", lw=1.2, label="Trench axis (PB2002)"),
        Line2D([], [], color=OI["blue"], lw=1.2, label="Fracture zone (GSFML)"),
        matplotlib.patches.Patch(
            fc=OI["verm"], alpha=0.35, ec=OI["verm"], label="LIP polygon"
        ),
    ]
    h2, l2 = ax.get_legend_handles_labels()
    ax.legend(
        handles=h + h2, loc="lower left", bbox_to_anchor=(0.44, 0.02), fontsize=8
    ).set_zorder(30)
    save(fig, "fabric_overview")


def fig_axis(t=1):
    S = SERIES[t]
    a, b = window(t)
    fig = plt.figure(figsize=(7.1, 5.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 1])
    ax = fig.add_subplot(gs[0])
    m = basemap(ax, lon0=160, lon1=216, lat0=47, lat1=62, res="i", dpar=5, dmer=10)
    x, y = m(S["axis_lon"], S["axis_lat"])
    xp, yp = m(S["pick_lon"], S["pick_lat"])
    m.plot(x, y, color="black", lw=1.2, label="PB2002 guide axis")
    m.plot(xp, yp, color=OI["verm"], lw=2.0, label="Thalweg pick (Viterbi)")
    import extract_series as es

    dn = np.column_stack([S["axis_lon"], S["axis_lat"]])
    n = len(dn)
    for i in range(0, n, 20):
        brg = es.bearing(dn[max(i - 2, 0)], dn[min(i + 2, n - 1)])
        p1 = es.dest(dn[i, 0], dn[i, 1], brg + np.pi / 2, es.UMAX)
        p2 = es.dest(dn[i, 0], dn[i, 1], brg - np.pi / 2, es.UMAX)
        xx, yy = m([p1[0], p2[0]], [p1[1], p2[1]])
        ax.plot(xx, yy, color="0.15", lw=0.6)
    ax.legend(loc="upper right", fontsize=8)
    tag(ax, "a")
    ax2 = fig.add_subplot(gs[1])
    s = S["s_km"]
    d = -S["d_m"] / 1000
    valid = (S["d_m"] <= an.DMAX) & (S["rise_land_m"] >= an.RLMIN)
    ax2.plot(s, d, color="0.55", lw=1.2, label="Raw thalweg depth")
    dd, _, _ = an.features(S, a, b, an.SMOOTH)
    ax2.plot(
        s[a : b + 1], -dd / 1000, color=OI["blue"], lw=2.0, label="Smoothed (15 km)"
    )
    ax2.scatter(
        s[~valid],
        d[~valid],
        s=10,
        marker="x",
        color=OI["verm"],
        lw=1.0,
        zorder=5,
        label="Trench-floor criterion failed",
    )
    ax2.set_xlim(s[0], s[-1])
    ax2.invert_yaxis()
    ax2.set_xlabel("Along-strike distance (km)")
    ax2.set_ylabel("Axial depth (km)")
    style(ax2, grid="y")
    ax2.legend(loc="upper center", fontsize=8, ncol=3)
    tag(ax2, "b", "left")
    save(fig, "axis_extraction")


def fig_cpmodel(t=1):
    S = SERIES[t]
    a, b = window(t)
    d, W, A = an.features(S, a, b, an.SMOOTH)
    s = S["s_km"][a : b + 1] - S["s_km"][a]
    pf = 2.0 ** np.arange(-2, 3.01, 0.5)
    counts = [len(an.segment(d, W, A, penf=p)[0]) for p in pf]
    cps, beta, Y = an.segment(d, W, A)
    fig = plt.figure(figsize=(7.1, 4.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 2.6])
    ax = fig.add_subplot(gs[0])
    ax.plot(pf, counts, "o-", color=OI["blue"], lw=1.2, ms=4)
    ax.axvline(1.0, color=OI["verm"], lw=0.7, ls="--")
    ax.annotate(
        f"BIC choice\nβ = {beta:.3f}\n{len(cps)} boundaries",
        xy=(1.0, len(cps)),
        xytext=(2.2, max(counts) * 0.82),
        fontsize=8.5,
        ha="center",
        arrowprops=dict(arrowstyle="-|>", lw=1.0),
    )
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.25, 0.5, 1, 2, 4, 8])
    ax.set_xticklabels(["¼", "½", "1", "2", "4", "8"])
    ax.set_xlabel("Penalty factor (× BIC β)")
    ax.set_ylabel("Number of boundaries")
    style(ax)
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    tag(ax, "a", "right")
    ax2 = fig.add_subplot(gs[1])
    off = [0, -4.5, -9]
    edges = [0] + cps + [len(d)]
    for k, (lab, col) in enumerate(
        (
            ("Axial depth d\n(up = shallower)", OI["blue"]),
            ("Width W", OI["orange"]),
            ("Asymmetry A", OI["purple"]),
        )
    ):
        ax2.plot(s, Y[:, k] + off[k], color="0.6", lw=1.2)
        for e0, e1 in zip(edges[:-1], edges[1:]):
            ax2.plot(
                s[e0:e1],
                np.full(e1 - e0, Y[e0:e1, k].mean()) + off[k],
                color=col,
                lw=2.0,
            )
        ax2.annotate(
            lab,
            xy=(1.015, off[k]),
            xycoords=("axes fraction", "data"),
            color=col,
            fontsize=8.5,
            va="center",
            ha="left",
            annotation_clip=False,
        )
    for c in cps:
        ax2.axvline(s[c] - 2.5, color=OI["verm"], lw=0.7, ls="--", alpha=0.8)
    ax2.set_xlim(0, s[-1])
    ax2.set_yticks([])
    ax2.set_xlabel("Along-strike distance (km)")
    ax2.set_ylabel("Standardised series (offset)")
    style(ax2, grid="x")
    ax2.tick_params(axis="y", which="both", left=False, right=False)
    tag(ax2, "b")
    save(fig, "changepoint_model")


def fig_panel():
    fig, axs = plt.subplots(5, 4, figsize=(7.1, 8.6), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)
    cat = {}
    for r in CAT:
        cat.setdefault(int(r["trench_id"]), []).append(
            (float(r["s_km"]), float(r["stability"]))
        )
    con = {}
    for r in CONTACTS:
        con.setdefault(int(r["trench_id"]), []).append(r)
    for t, ax in zip(range(1, 21), axs.flat):
        S = SERIES[t]
        s = S["s_km"]
        d = -S["d_m"] / 1000
        w = window(t)
        ax.plot(s, d, color="0.6", lw=1.0)
        if w:
            a, b = w
            ax.axvspan(s[a], s[b], color=OI["sky"], alpha=0.12, lw=0)
            dd, _, _ = an.features(S, a, b, an.SMOOTH)
            ax.plot(s[a : b + 1], -dd / 1000, color=OI["blue"], lw=1.2)
            for x, st in cat.get(t, []):
                ax.axvline(
                    x,
                    color=OI["verm"],
                    lw=0.7 if st >= an.STAB_MIN else 0.5,
                    ls="-" if st >= an.STAB_MIN else ":",
                    alpha=0.9,
                )
            yb = d.min() - 0.06 * (d.max() - d.min() + 0.5)
            for r in con.get(t, []):
                lab, col, mk = CLS[r["cls"]]
                if r["cls"] == "plateau":
                    ax.plot(
                        [float(r["s1_km"]), float(r["s2_km"])],
                        [yb, yb],
                        color=col,
                        lw=3,
                        solid_capstyle="butt",
                    )
                else:
                    ax.plot(
                        float(r["s1_km"]),
                        yb,
                        marker=mk,
                        color=col,
                        ms=4,
                        mec="black",
                        mew=0.4,
                        ls="none",
                    )
            ttl = f"{t} {S['name']}"
        else:
            ttl = f"{t} {S['name']} (excluded)"
        ax.set_title(ttl, fontsize=8.5, loc="left", pad=2)
        ax.set_xlim(s[0], s[-1])
        style(ax, grid="y")
        rng_ = d.max() - d.min() + 0.5
        ax.set_ylim(d.max() + 0.03 * rng_, d.min() - 0.12 * rng_)
        ax.tick_params(labelsize=8)
    for ax in axs[-1]:
        ax.set_xlabel("Distance (km)", fontsize=8.5)
    for ax in axs[:, 0]:
        ax.set_ylabel("Depth (km)", fontsize=8.5)
    h = [
        Line2D([], [], color=OI["blue"], lw=1.2, label="Smoothed axial depth"),
        Line2D([], [], color=OI["verm"], lw=0.7, label="Stable boundary"),
        Line2D([], [], color=OI["verm"], lw=0.5, ls=":", label="Other boundary"),
    ]
    h += [
        Line2D(
            [],
            [],
            color=CLS["fracture_zone"][1],
            marker="s",
            ls="none",
            mec="black",
            mew=0.4,
            label="Fracture zone contact",
        ),
        Line2D(
            [],
            [],
            color=CLS["seamount_ridge"][1],
            marker="^",
            ls="none",
            mec="black",
            mew=0.4,
            label="Seamount contact",
        ),
        Line2D(
            [],
            [],
            color=CLS["plateau"][1],
            lw=3,
            label="LIP plateau / ridge entry interval",
        ),
    ]
    fig.legend(
        handles=h,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=3,
        fontsize=8.5,
        frameon=False,
    )
    save(fig, "alongstrike_panel")


def fig_boundary_map():
    fig, ax = plt.subplots(figsize=(7.1, 4.6))
    m = basemap(ax, relief="grey")
    for t in range(1, 21):
        S = SERIES[t]
        window(t)
        x, y = m(S["pick_lon"], S["pick_lat"])
        m.plot(x, y, color="0.2", lw=0.8, zorder=3)
    con = {}
    for r in CONTACTS:
        con.setdefault(int(r["trench_id"]), {}).setdefault(r["cls"], []).append(
            (float(r["s1_rel_km"]), float(r["s2_rel_km"]))
        )
    xs, ys, cls_hit = [], [], []
    for r in CAT:
        t = int(r["trench_id"])
        x, y = m(float(r["lon"]), float(r["lat"]))
        hit = None
        for c in ("plateau", "fracture_zone", "seamount_ridge"):
            iv = con.get(t, {}).get(c)
            if iv and an.dist_to([float(r["s_rel_km"])], iv)[0] <= an.D_TOL:
                hit = c
                break
        xs.append(x)
        ys.append(y)
        cls_hit.append(hit)
    xs = np.array(xs)
    ys = np.array(ys)
    nohit = np.array([h is None for h in cls_hit])
    ax.scatter(
        xs[nohit],
        ys[nohit],
        s=6,
        c="black",
        zorder=5,
        label=f"Other boundaries ({nohit.sum()})",
    )
    for c, (lab, col, mk) in CLS.items():
        sel = np.array([h == c for h in cls_hit])
        if sel.any():
            ax.scatter(
                xs[sel],
                ys[sel],
                s=34,
                marker=mk,
                c=col,
                edgecolors="black",
                linewidths=0.5,
                zorder=7,
                label=f"≤15 km of {HITLAB[c]} ({sel.sum()})",
            )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.055),
        fontsize=8,
        ncol=2,
        columnspacing=1.5,
        handletextpad=0.3,
        frameon=False,
    )
    save(fig, "boundary_map")


def fig_coincidence():
    nulls = json.load(open(os.path.join(HERE, "results", "null_distributions.json")))
    rows = {(r["set"], r["cls"]): r for r in COIN}
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.6), constrained_layout=True)
    for k, (c, ax) in enumerate(zip(CLS, axs)):
        nd = nulls[f"all_{c}"]
        r = rows[("all", c)]
        lab, col, mk = CLS[c]
        u, cs = np.array(nd["uni"]), np.array(nd["cs"])
        lo, hi = min(u.min(), cs.min(), nd["obs"]), max(u.max(), cs.max(), nd["obs"])
        bins = np.linspace(lo, hi, 40)
        ax.hist(u, bins=bins, color="0.7", alpha=0.8, label="Uniform null")
        ax.hist(
            cs,
            bins=bins,
            histtype="step",
            color=col,
            lw=1.2,
            label="Circular-shift null (class colour)",
        )
        ax.axvline(nd["obs"], color="black", lw=2.0, label="Observed")
        fp = lambda p: "> 0.99" if p > 0.99 else f"= {p:.2f}"
        ax.set_title(
            f"{lab}\np {fp(float(r['p_uni']))} (uniform), p {fp(float(r['p_cs']))} (shift)",
            fontsize=8.5,
        )
        ax.set_xlabel("Mean nearest distance (km)")
        if k == 0:
            ax.set_ylabel("Count")
        style(ax, grid="y")
        tag(ax, "abc"[k])
    fig.legend(
        *axs[0].get_legend_handles_labels(),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=3,
        fontsize=8.5,
        frameon=False,
    )
    save(fig, "coincidence")


def fig_classfrac():
    rows = {(r["set"], r["cls"]): r for r in COIN}
    fig, axs = plt.subplots(
        1, 3, figsize=(7.1, 2.5), constrained_layout=True, sharey=False
    )
    for k, (c, ax) in enumerate(zip(CLS, axs)):
        lab, col, mk = CLS[c]
        for j, st in enumerate(("all", "stable")):
            r = rows.get((st, c))
            if r is None:
                continue
            tol = np.array(an.TOLS) + (j - 0.5) * 4
            ob = [float(r[f"frac{int(t)}"]) for t in an.TOLS]
            nu = [float(r[f"frac{int(t)}_uni"]) for t in an.TOLS]
            hi = [float(r[f"frac{int(t)}_uni_hi"]) for t in an.TOLS]
            lo = [float(r[f"frac{int(t)}_uni_lo"]) for t in an.TOLS]
            ax.vlines(tol, lo, hi, color="0.5", lw=1.2)
            ax.plot(tol, nu, "_", color="0.3", ms=8, mew=1.2)
            ax.plot(
                tol,
                ob,
                marker=mk if st == "all" else "o",
                ls="none",
                color=col,
                mfc=col if st == "all" else "white",
                mec=col,
                ms=6,
                mew=1.2,
                label="All boundaries" if st == "all" else "Stable boundaries",
            )
        ax.set_xticks(an.TOLS)
        ax.set_xlabel("Tolerance (km)")
        if k == 0:
            ax.set_ylabel("Fraction of boundaries")
        ax.set_title(lab, fontsize=9)
        style(ax, grid="y")
        tag(ax, "abc"[k])
        ax.set_ylim(bottom=0)
    hh = [
        Line2D(
            [],
            [],
            marker="s",
            ls="none",
            color="0.3",
            ms=6,
            label="All boundaries (filled)",
        ),
        Line2D(
            [],
            [],
            marker="o",
            ls="none",
            mfc="white",
            mec="0.3",
            ms=6,
            mew=1.2,
            label="Stable boundaries (open)",
        ),
        Line2D(
            [],
            [],
            marker="_",
            color="0.3",
            ms=8,
            mew=1.2,
            ls="none",
            label="Uniform-null mean, 95% range (bar)",
        ),
    ]
    fig.legend(
        handles=hh,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=3,
        fontsize=8.5,
        frameon=False,
    )
    save(fig, "classfrac")


def fig_contrast():
    rng = np.random.default_rng(an.SEED)
    obs = {k: [] for k in "dWA"}
    ref = {k: [] for k in "dWA"}
    for r in CAT:
        for k in "dWA":
            obs[k].append(abs(float(r[f"{k}z_right"]) - float(r[f"{k}z_left"])))
    for t in range(1, 21):
        w = window(t)
        if not w:
            continue
        a, b = w
        d, W, A = an.features(SERIES[t], a, b, an.SMOOTH)
        Y = np.column_stack([an.zs(d), an.zs(W), an.zs(A)])
        n = len(Y)
        nb = sum(1 for r in CAT if int(r["trench_id"]) == t)
        for _ in range(50):
            pos = np.sort(
                rng.choice(np.arange(an.MINSIZE, n - an.MINSIZE), nb, replace=False)
            )
            e = np.r_[0, pos, n]
            for i, c in enumerate(pos):
                l, rr = Y[e[i] : c], Y[c : e[i + 2]]
                for j, k in enumerate("dWA"):
                    ref[k].append(abs(rr[:, j].mean() - l[:, j].mean()))
    fig, axs = plt.subplots(
        1, 3, figsize=(7.1, 2.5), constrained_layout=True, sharey=True
    )
    names = {
        "d": ("Axial depth", OI["blue"]),
        "W": ("Width", OI["orange"]),
        "A": ("Asymmetry", OI["purple"]),
    }
    bins = np.linspace(0, 3, 31)
    for k, (key, ax) in enumerate(zip("dWA", axs)):
        lab, col = names[key]
        ax.hist(
            ref[key], bins=bins, density=True, color="0.75", label="Random positions"
        )
        ax.hist(
            obs[key],
            bins=bins,
            density=True,
            histtype="step",
            color=col,
            lw=2.0,
            label="Detected boundaries",
        )
        mo, mr = np.median(obs[key]), np.median(ref[key])
        ax.set_title(f"{lab}: median {mo:.2f} vs {mr:.2f} σ", fontsize=8.5)
        ax.set_xlabel("|Step| (standard deviations)")
        style(ax, grid="y")
        tag(ax, "abc"[k], "right")
        if k == 0:
            ax.set_ylabel("Density")
        ax.set_xlim(0, 3)
    axs[1].legend(loc="center right", fontsize=8)
    save(fig, "segment_contrast")
    return {
        k: (
            float(np.median(obs[k])),
            float(np.median(ref[k])),
            float(np.mean(np.array(obs[k]) > 0.5)),
        )
        for k in "dWA"
    }


def fig_sensitivity(t=1):
    rows = [r for r in rcsv("sensitivity_grid.csv") if int(r["trench_id"]) == t]
    fig = plt.figure(figsize=(7.1, 4.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 2.2], height_ratios=[1, 1])
    ax = fig.add_subplot(gs[:, 0])
    cols = {1: OI["sky"], 3: OI["blue"], 5: OI["black"]}
    mks = {1: "o", 3: "s", 5: "^"}
    for sm in an.GRID_SMOOTH:
        for md, ls in (("l2", "-"), ("normal", "--")):
            rr = [
                r
                for r in rows
                if int(r["smooth"]) == sm
                and r["cost"] == md
                and int(r["min_size"]) == an.MINSIZE
            ]
            rr.sort(key=lambda r: float(r["pen_factor"]))
            ax.plot(
                [float(r["pen_factor"]) for r in rr],
                [int(r["n_boundaries"]) for r in rr],
                ls=ls,
                marker=mks[sm],
                color=cols[sm],
                lw=1.2,
                ms=4,
                label=f"{str(sm * 5) + ' km' if sm > 1 else 'none'}, {'mean' if md == 'l2' else 'mean+cov'}",
            )
    ax.set_xscale("log", base=2)
    ax.set_xticks(an.GRID_PENF)
    ax.set_xticklabels(["0.5", "1", "2"])
    ax.set_xlabel("Penalty factor (× BIC)")
    ax.set_ylabel("Number of boundaries")
    style(ax)
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    fig.legend(
        *ax.get_legend_handles_labels(),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=3,
        fontsize=8.5,
        frameon=False,
        title="Smoothing, cost function",
        title_fontsize=8.5,
    )
    tag(ax, "a", "right")
    S = SERIES[t]
    a, b = window(t)
    s = S["s_km"]
    cat = [r for r in CAT if int(r["trench_id"]) == t]
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1], sharex=ax2)
    ax2.plot(s[a : b + 1], -S["d_m"][a : b + 1] / 1000, color="0.5", lw=1.2)
    ax2.invert_yaxis()
    ax2.set_ylabel("Depth (km)")
    style(ax2, grid="y")
    tag(ax2, "b")
    x = np.array([float(r["s_km"]) for r in cat])
    st = np.array([float(r["stability"]) for r in cat])
    ax3.bar(x, st, width=18, color=np.where(st >= an.STAB_MIN, OI["verm"], "0.6"))
    ax3.axhline(an.STAB_MIN, color="black", lw=0.7, ls="--")
    ax3.set_ylim(0, 1.05)
    ax3.set_xlim(s[a], s[b])
    ax3.set_xlabel("Along-strike distance (km)")
    ax3.set_ylabel("Recurrence fraction")
    style(ax3, grid="y")
    tag(ax3, "c")
    save(fig, "sensitivity")


if __name__ == "__main__":
    import matplotlib.patheffects
    import matplotlib.patches

    fig_studyarea()
    fig_fabric()
    fig_axis()
    fig_cpmodel()
    fig_panel()
    fig_boundary_map()
    fig_coincidence()
    fig_classfrac()
    print("contrast", fig_contrast())
    fig_sensitivity()
