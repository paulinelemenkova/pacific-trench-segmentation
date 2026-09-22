#!/usr/bin/env python3
import math
import numpy as np
import pandas as pd
import pygmt

AXIS_FILE = "axes_full/1_Aleutian.txt"
GRID = "@earth_relief_02m"
DS, L, DU, SMOOTH = 12.0, 50.0, 4.0, 9
MAP_REGION = [162, 217, 49.5, 61]
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

def build_geometry():
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
    offs = np.array([-L + k * DU for k in range(int(2 * L // DU) + 1)])
    lons, lats = [], []
    for i, (lon, lat) in enumerate(dense):
        for u in offs:
            b = brg[i] + math.pi / 2 if u >= 0 else brg[i] - math.pi / 2
            plon, plat = dest(lon, lat, b, abs(u))
            lons.append(plon); lats.append(plat)
    return np.array(dense), np.array(lons), np.array(lats), offs

def main():
    dense, lons, lats, offs = build_geometry()
    n, no = len(dense), len(offs)

    sampled = pygmt.grdtrack(points=pd.DataFrame({"x": lons, "y": lats}),
                             grid=GRID, newcolname="z")
    z = sampled["z"].to_numpy().reshape(n, no)
    plon = lons.reshape(n, no); plat = lats.reshape(n, no)

    kmin = np.nanargmin(z, axis=1)
    d = z[np.arange(n), kmin] / 1000.0
    picks = np.column_stack([plon[np.arange(n), kmin], plat[np.arange(n), kmin]])
    s = np.arange(n) * DS
    dsm = np.convolve(d, np.ones(SMOOTH) / SMOOTH, mode="same")

    np.savetxt("axis_dense.txt", dense, fmt="%.4f")
    np.savetxt("picks.txt", picks, fmt="%.4f")
    np.savetxt("series_raw.txt", np.column_stack([s, d]), fmt="%.2f")
    np.savetxt("series_smooth.txt", np.column_stack([s, dsm]), fmt="%.2f")
    step = max(1, n // 26)
    with open("fan.txt", "w") as f:
        for i in range(0, n, step):
            f.write(f"> \n{plon[i,0]} {plat[i,0]}\n{plon[i,-1]} {plat[i,-1]}\n")
    open("legend.txt", "w").write(
        "N 3\n"
        "S 0.4c - 0.9c - 1.4p,gold1,-- 1.1c Digitised axis\n"
        "S 0.4c - 0.9c - 1.8p,red 1.1c Deepest-point locus (Eq. 1)\n"
        "S 0.4c - 0.9c - 0.6p,chartreuse1 1.1c Cross-profiles (fan)\n")

    smax = float(s[-1]); zlo = round(float(d.min()) - 0.3, 2); zhi = round(float(d.max()) + 0.3, 2)

    fig = pygmt.Figure()
    pygmt.config(FONT_TITLE="11p,Helvetica", FONT_ANNOT_PRIMARY="9p",
                 FONT_LABEL="10p", MAP_FRAME_TYPE="plain",
                 MAP_GRID_PEN_PRIMARY="0.25p,gray80")

    fig.basemap(region=[0, smax, zlo, zhi], projection="X18c/5c",
                frame=["WSne+t(b) Along-strike axial-depth series d(s)",
                       "xa500f100g500+lAlong-strike distance (km)",
                       "ya1f0.5g1+lAxial depth (km)"])
    fig.plot(data="series_raw.txt", pen="0.6p,gray70")
    fig.plot(data="series_smooth.txt", pen="1p,navy")

    fig.shift_origin(yshift="7.0c")
    pygmt.makecpt(cmap="geo", series=[-8000, 5000, 100])
    with pygmt.config(MAP_GRID_PEN_PRIMARY="0.2p,white"):
        fig.grdimage(grid=GRID, region=MAP_REGION, projection="M18c", shading="+d",
                     frame=["xa10f5g10", "ya4f2g4",
                            "WsNe+t(a) Cross-profile sampling and deepest-point picks"])
    fig.coast(land="gray82", shorelines="0.3p,gray45", resolution="i")
    fig.plot(data="fan.txt", pen="0.5p,chartreuse1")
    fig.plot(data="axis_dense.txt", pen="1.4p,gold1,--")
    fig.plot(data="picks.txt", pen="1.8p,red")
    fig.legend(spec="legend.txt", position="JBC+jTC+w15c+o0/0.35c",
               box="+gwhite@10+p0.5p,gray50")

    fig.savefig("fig_axis_extraction.pdf")
    fig.savefig("fig_axis_extraction.png", dpi=300)
    print("wrote fig_axis_extraction.pdf and fig_axis_extraction.png")

if __name__ == "__main__":
    main()
