#!/usr/bin/env python3
import math
import statistics

import pandas as pd
import pygmt

AXIS_FILE = "axes_full/1_Aleutian.txt"
GRID = "@earth_relief_02m"
DS, L, DU, LMIN = 15.0, 50.0, 5.0, 4
SMOOTHS = [3, 5, 7, 9]
FACS = [1, 1.5, 2, 2.7, 3.5, 4.5, 6, 7.5, 9, 11]
THR, BINW = 0.6, 30.0
COLORS = {3: "gray55", 5: "dodgerblue", 7: "navy", 9: "firebrick"}
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

def mv(a, w):
    h = w // 2
    return [sum(a[max(0, i - h):min(len(a), i + h + 1)]) /
            len(a[max(0, i - h):min(len(a), i + h + 1)]) for i in range(len(a))]

def partition(dd, beta):
    n = len(dd); c1 = [0.0] * (n + 1); c2 = [0.0] * (n + 1)
    for k in range(n):
        c1[k + 1] = c1[k] + dd[k]; c2[k + 1] = c2[k] + dd[k] ** 2

    def sc(a, b):
        m2 = b - a; s1 = c1[b] - c1[a]; s2 = c2[b] - c2[a]
        return s2 - s1 * s1 / m2

    F = [None] * (n + 1); F[0] = -beta; prev = [0] * (n + 1)
    for t in range(1, n + 1):
        hi = t - LMIN
        cs = range(0, hi + 1) if hi >= 0 else [0]
        bv = None; bs = 0
        for a in cs:
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

def axial_series():
    ax = [tuple(map(float, ln.split()[:2])) for ln in open(AXIS_FILE)
          if ln.strip() and not ln.startswith("#")]
    cum = [0.0]
    for i in range(1, len(ax)):
        cum.append(cum[-1] + hav(ax[i - 1], ax[i]))
    dense = []
    for s in [k * DS for k in range(int(cum[-1] // DS) + 1)]:
        j = 0
        while j < len(cum) - 2 and cum[j + 1] < s:
            j += 1
        t = (s - cum[j]) / max(cum[j + 1] - cum[j], 1e-9)
        dense.append((ax[j][0] + t * (ax[j + 1][0] - ax[j][0]),
                      ax[j][1] + t * (ax[j + 1][1] - ax[j][1])))
    n = len(dense)
    brg = [bearing(dense[max(i - 1, 0)], dense[min(i + 1, n - 1)]) for i in range(n)]
    offs = [-L + k * DU for k in range(int(2 * L // DU) + 1)]
    lons, lats, idx = [], [], []
    for i, (lon, lat) in enumerate(dense):
        for u in offs:
            b = brg[i] + math.pi / 2 if u >= 0 else brg[i] - math.pi / 2
            plon, plat = dest(lon, lat, b, abs(u))
            lons.append(plon); lats.append(plat); idx.append(i)
    sampled = pygmt.grdtrack(points=pd.DataFrame({"x": lons, "y": lats}),
                             grid=GRID, newcolname="z")
    zc = sampled["z"].to_list()
    best = {}
    for i, z in zip(idx, zc):
        if z == z and (i not in best or z < best[i]):
            best[i] = z
    s = [i * DS for i in range(n) if i in best]
    d = [best[i] / 1000.0 for i in range(n) if i in best]
    return s, d

def main():
    s, d = axial_series()
    n = len(d); smax = s[-1]

    allbounds = []
    for w in SMOOTHS:
        dsm = mv(d, w); mu = sum(dsm) / n
        sd = (sum((x - mu) ** 2 for x in dsm) / n) ** 0.5 or 1.0
        dd = [(x - mu) / sd for x in dsm]
        sig2 = 0.5 * statistics.pvariance([dd[k] - dd[k - 1] for k in range(1, n)])
        pts = []
        for f in FACS:
            beta = f * max(sig2, 1e-6) * math.log(n)
            cps = partition(dd, beta); pts.append((f, len(cps)))
            allbounds += [s[c] for c in cps]
        open(f"lineK_{w}.txt", "w").write("\n".join(f"{f} {K}" for f, K in pts) + "\n")

    ncomb = len(SMOOTHS) * len(FACS); nb = int(smax // BINW) + 1
    cnt = [0] * nb
    for sb in allbounds:
        cnt[min(int(sb // BINW), nb - 1)] += 1
    rec = [c / ncomb for c in cnt]
    open("recur.txt", "w").write("\n".join(f"{(k+0.5)*BINW} {rec[k]}" for k in range(nb)) + "\n")
    stable = []; k = 0
    while k < nb:
        if rec[k] >= THR:
            j = k
            while j + 1 < nb and rec[j + 1] >= THR:
                j += 1
            bb = max(range(k, j + 1), key=lambda b: rec[b]); stable.append((bb + 0.5) * BINW); k = j + 1
        else:
            k += 1
    open("stable.txt", "w").write("\n".join(f"{x} 1.03" for x in stable) + "\n")
    d7 = mv(d, 7)
    open("profileB.txt", "w").write("\n".join(f"{s[k]} {d7[k]}" for k in range(n)) + "\n")
    with open("legA.txt", "w") as f:
        f.write("L 9p,Helvetica C Smoothing width (nodes)\nN 4\n")
        for w in SMOOTHS:
            f.write(f"S 0.3c - 1.0c - 1p,{COLORS[w]} 0.4c {w}\n")
    dmin, dmax = min(d) - 0.3, max(d) + 0.3
    kmax = max(K for w in SMOOTHS for _, K in
               [tuple(map(float, ln.split())) for ln in open(f"lineK_{w}.txt")]) + 2

    fig = pygmt.Figure()
    pygmt.config(FONT_TITLE="11p,Helvetica", FONT_ANNOT_PRIMARY="9p", FONT_LABEL="10p",
                 MAP_FRAME_TYPE="plain", MAP_GRID_PEN_PRIMARY="0.2p,gray88")

    fig.basemap(region=[0, smax, 0, 1.1], projection="X17c/5c",
                frame=["WSne+t(b) Boundary recurrence across the parameter grid (Aleutian)",
                       "xa500f100g500+lAlong-strike distance (km)",
                       "ya0.2f0.1g0.2+lRecurrence (fraction of grid)"])
    fig.plot(data="recur.txt", style="b25u+b0", fill="steelblue", pen="0.2p,gray30")
    fig.plot(x=[0, smax], y=[THR, THR], pen="0.4p,red,--")
    fig.plot(data="stable.txt", style="i0.30c", fill="red", pen="0.3p,black", no_clip=True)
    fig.text(x=1300, y=THR, text=f"stable (>={THR})", font="8p,Helvetica-Bold,red",
             justify="CB", offset="0/0.12c", fill="white@25", no_clip=True)
    fig.basemap(region=[0, smax, dmin, dmax], projection="X17c/5c",
                frame=["E", "ya2f1+lAxial depth (km)"])
    fig.plot(data="profileB.txt", pen="0.8p,gray45")

    fig.shift_origin(yshift="7.3c")
    fig.basemap(region=[FACS[0], FACS[-1], 0, kmax], projection="X17c/4.2c",
                frame=["WSne+t(a) Number of boundaries vs penalty, by smoothing width",
                       "xa2f1g1+lPenalty factor (@~b@~ / BIC)", "ya5f1g5+lBoundaries K"])
    for w in SMOOTHS:
        fig.plot(data=f"lineK_{w}.txt", pen=f"0.5p,{COLORS[w]}")
        fig.plot(data=f"lineK_{w}.txt", style="c0.08c", fill=COLORS[w])
    fig.legend(spec="legA.txt", position="jBL+w9c+o0.3c/0.3c",
               box="+gwhite@10+p0.4p,gray50")

    fig.savefig("fig_sensitivity.pdf")
    fig.savefig("fig_sensitivity.png", dpi=300)
    print("wrote fig_sensitivity.pdf and fig_sensitivity.png")

if __name__ == "__main__":
    main()
