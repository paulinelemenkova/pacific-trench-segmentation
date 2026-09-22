#!/usr/bin/env python3
import glob
import math
import os
import re
import statistics

import pandas as pd
import pygmt

GRID = "@earth_relief_04m"
DS, L, DU = 15.0, 40.0, 5.0
SMOOTH, PEN_FAC, LMIN = 7, 4.5, 5
R = 6371.0

def hav(a, b):
    lo1, la1 = map(math.radians, a); lo2, la2 = map(math.radians, b)
    h = math.sin((la2 - la1) / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin((lo2 - lo1) / 2) ** 2
    return 2 * R * math.asin(math.sqrt(h))

def bearing(a, b):
    lo1, la1 = map(math.radians, a); lo2, la2 = map(math.radians, b)
    y = math.sin(lo2 - lo1) * math.cos(la2)
    x = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(lo2 - lo1)
    return math.atan2(y, x)

def dest(lon, lat, brg, dist):
    la1, lo1, dR = math.radians(lat), math.radians(lon), dist / R
    la2 = math.asin(math.sin(la1) * math.cos(dR) + math.cos(la1) * math.sin(dR) * math.cos(brg))
    lo2 = lo1 + math.atan2(math.sin(brg) * math.sin(dR) * math.cos(la1),
                           math.cos(dR) - math.sin(la1) * math.sin(la2))
    return (math.degrees(lo2) % 360, math.degrees(la2))

def movavg(a, w):
    h = w // 2
    return [sum(a[max(0, i - h):min(len(a), i + h + 1)]) /
            len(a[max(0, i - h):min(len(a), i + h + 1)]) for i in range(len(a))]

def partition(dd, beta):
    n = len(dd); c1 = [0.0] * (n + 1); c2 = [0.0] * (n + 1)
    for k in range(n):
        c1[k + 1] = c1[k] + dd[k]; c2[k + 1] = c2[k] + dd[k] ** 2

    def sc(a, b):
        nn = b - a; s1 = c1[b] - c1[a]; s2 = c2[b] - c2[a]
        return s2 - s1 * s1 / nn

    F = [None] * (n + 1); F[0] = -beta; prev = [0] * (n + 1)
    for t in range(1, n + 1):
        hi = t - LMIN
        cands = range(0, hi + 1) if hi >= 0 else [0]
        bv = None; bs = 0
        for a in cands:
            cc = F[a] + sc(a, t) + beta
            if bv is None or cc < bv:
                bv = cc; bs = a
        F[t] = bv; prev[t] = bs
    cps = []; t = n
    while t > 0:
        a = prev[t]
        if a > 0:
            cps.append(a)
        t = a
    return sorted(cps)

def load_trenches():
    TR = {}
    for p in glob.glob("axes_full/*.txt"):
        m = re.match(r"(\d+)_(.+)\.txt$", os.path.basename(p))
        if not m:
            continue
        TR[int(m.group(1))] = (m.group(2).replace("_", " "),
            [tuple(map(float, ln.split()[:2])) for ln in open(p)
             if ln.strip() and not ln.startswith("#")])
    return TR

def main():
    TR = load_trenches()
    lons, lats, tids, nids = [], [], [], []
    for tid in sorted(TR):
        _, AX = TR[tid]
        cum = [0.0]
        for i in range(1, len(AX)):
            cum.append(cum[-1] + hav(AX[i - 1], AX[i]))
        dense = []
        for s in [k * DS for k in range(int(cum[-1] // DS) + 1)]:
            j = 0
            while j < len(cum) - 2 and cum[j + 1] < s:
                j += 1
            t = (s - cum[j]) / max(cum[j + 1] - cum[j], 1e-9)
            dense.append((AX[j][0] + t * (AX[j + 1][0] - AX[j][0]),
                          AX[j][1] + t * (AX[j + 1][1] - AX[j][1])))
        n = len(dense)
        brg = [bearing(dense[max(i - 1, 0)], dense[min(i + 1, n - 1)]) for i in range(n)]
        offs = [-L + k * DU for k in range(int(2 * L // DU) + 1)]
        for i, (lon, lat) in enumerate(dense):
            for u in offs:
                b = brg[i] + math.pi / 2 if u >= 0 else brg[i] - math.pi / 2
                plon, plat = dest(lon, lat, b, abs(u))
                lons.append(plon); lats.append(plat); tids.append(tid); nids.append(i)

    sampled = pygmt.grdtrack(points=pd.DataFrame({"x": lons, "y": lats}),
                             grid=GRID, newcolname="z")
    zc = sampled["z"].to_list()
    best = {}
    for tid, nid, z in zip(tids, nids, zc):
        if z == z:
            key = (tid, nid)
            if key not in best or z < best[key]:
                best[key] = z

    panels = []
    for tid in sorted(TR):
        name = TR[tid][0]
        ns = sorted(k[1] for k in best if k[0] == tid)
        if not ns:
            continue
        s = [nid * DS for nid in ns]; d = [best[(tid, nid)] / 1000.0 for nid in ns]
        d = movavg(d, SMOOTH); n = len(d)
        mu = sum(d) / n; sd = (sum((x - mu) ** 2 for x in d) / n) ** 0.5 or 1.0
        dd = [(x - mu) / sd for x in d]
        sig2 = 0.5 * statistics.pvariance([dd[k] - dd[k - 1] for k in range(1, n)]) if n > 2 else 1.0
        beta = PEN_FAC * max(sig2, 1e-6) * math.log(max(n, 2))
        cps = partition(dd, beta)
        open(f"series_{tid}.txt", "w").write("\n".join(f"{s[k]} {d[k]}" for k in range(n)) + "\n")
        dmin = min(d) - 0.3; dmax = max(d) + 0.3
        open(f"bounds_{tid}.txt", "w").write(
            "\n".join(f"> \n{s[c]} {dmin}\n{s[c]} {dmax}" for c in cps) + "\n")
        panels.append((tid, name, max(s[-1], 1), dmin, dmax))
    open("legend.txt", "w").write(
        "N 2\n"
        "S 0.7c - 1.1c - 0.8p,navy 1.4c Axial-depth profile d(s)\n"
        "S 0.7c - 1.1c - 0.4p,red,-- 1.4c Detected segment boundary\n")

    fig = pygmt.Figure()
    pygmt.config(FONT_HEADING="12p,Helvetica", FONT_TITLE="9p,Helvetica",
                 FONT_ANNOT_PRIMARY="6p", FONT_LABEL="7p", MAP_FRAME_TYPE="plain",
                 MAP_FRAME_PEN="0.6p", MAP_TICK_LENGTH_PRIMARY="0.08c",
                 MAP_GRID_PEN_PRIMARY="0.15p,gray90")
    with fig.subplot(nrows=5, ncols=4, figsize=("17.85c", "17.65c"),
                     margins=["0.35c", "0.85c"],
                     title="Along-strike axial-depth profiles with detected segment boundaries"):
        for idx, (tid, name, smax, dmin, dmax) in enumerate(panels):
            with fig.set_panel(idx):
                fig.basemap(region=[0, smax, dmin, dmax], projection="X4.2c/2.85c",
                            frame=["xafg", "yafg", f"WSne+t{tid}. {name}"])
                fig.plot(data=f"bounds_{tid}.txt", pen="0.4p,red,--")
                fig.plot(data=f"series_{tid}.txt", pen="0.8p,navy")

    fig.legend(spec="legend.txt", position="JBC+jTC+w12c+o0/0.7c",
               box="+gwhite+p0.5p,gray50")

    fig.savefig("fig_alongstrike_panel.pdf")
    fig.savefig("fig_alongstrike_panel.png", dpi=300)
    print("wrote fig_alongstrike_panel.pdf and fig_alongstrike_panel.png")

if __name__ == "__main__":
    main()
