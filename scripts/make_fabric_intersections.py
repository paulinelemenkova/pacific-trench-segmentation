#!/usr/bin/env python3
import glob
import math
import os
import re
import csv

import numpy as np

CONFIG = {
    "fracture_zone": [("fabric/FZ.xy", "line")],
    "ridge":         [("fabric/seamounts.txt", "seamount")],

    "plateau":       [("fabric/lips.gmt", "polygon")],
}

MIN_SEAMOUNT_HEIGHT_M = 1000.0

BUFFER_KM = {"line": 50.0, "polygon": 90.0, "seamount": 60.0}
AXES_DIR  = "axes_full"
OUT_CSV   = "fabric_intersections.csv"
DEDUP_KM  = 25.0
AXIS_DS   = 5.0
FEAT_DS   = 5.0
R_EARTH   = 6371.0

def hav(lon1, lat1, lon2, lat2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi, dl = math.radians(lat2 - lat1), math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R_EARTH * math.asin(math.sqrt(a))

def densify(pts, ds):

    if len(pts) < 2:
        return list(pts)
    out = [pts[0]]
    for (lo1, la1), (lo2, la2) in zip(pts[:-1], pts[1:]):
        d = hav(lo1, la1, lo2, la2)
        k = max(int(d // ds), 1)
        for m in range(1, k + 1):
            t = m / k
            out.append((lo1 + t * (lo2 - lo1), la1 + t * (la2 - la1)))
    return out

def cum_dist(pts):
    s = [0.0]
    for (lo1, la1), (lo2, la2) in zip(pts[:-1], pts[1:]):
        s.append(s[-1] + hav(lo1, la1, lo2, la2))
    return s

def norm360(lon):
    return lon % 360.0

def _rows_from_csv(path):
    with open(path) as f:
        rdr = csv.reader(f)
        rows = list(rdr)
    if not rows:
        return []
    hdr = [c.strip().lower() for c in rows[0]]
    ix = iy = None
    for cand in ("lon", "longitude", "x"):
        if cand in hdr:
            ix = hdr.index(cand); break
    for cand in ("lat", "latitude", "y"):
        if cand in hdr:
            iy = hdr.index(cand); break
    start = 1 if (ix is not None and iy is not None) else 0
    if ix is None or iy is None:
        ix, iy = 0, 1
    out = []
    for r in rows[start:]:
        try:
            out.append((float(r[ix]), float(r[iy])))
        except (ValueError, IndexError):
            pass
    return [out]

def read_features(path, geom):

    if not os.path.exists(path):
        return None
    if path.lower().endswith(".csv"):
        pts = _rows_from_csv(path)[0]
        pts = [(norm360(lo), la) for lo, la in pts]
        return [[p] for p in pts] if geom == "point" else [pts]

    feats, cur = [], []
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        if ln.startswith(">"):
            if cur:
                feats.append(cur); cur = []
            continue
        parts = ln.replace(",", " ").split()
        try:
            lo, la = float(parts[0]), float(parts[1])
        except (ValueError, IndexError):
            continue
        cur.append((norm360(lo), la))
    if cur:
        feats.append(cur)
    if geom == "point":
        return [[p] for f in feats for p in f]
    return feats

def read_seamounts(path, min_h):

    if not os.path.exists(path):
        return None
    pts = []
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        p = ln.replace(",", " ").split()
        try:
            lo, la, h = float(p[0]), float(p[1]), float(p[5])
        except (ValueError, IndexError):
            continue
        if h >= min_h:
            pts.append((norm360(lo), la))
    return pts

def seamount_contacts(pts, axes, buffer):

    P = np.asarray(pts); N = len(P)
    best_d = np.full(N, np.inf); best_t = np.full(N, -1, int); best_k = np.zeros(N, int)
    for tid, axis in axes.items():
        ax = axis["pts"]
        kx = 111.320 * math.cos(math.radians(float(ax[:, 1].mean()))); ky = 110.574
        for a in range(0, N, 4000):
            Pb = P[a:a + len(P[a:a + 4000])]
            dx = (Pb[:, 0][:, None] - ax[:, 0][None, :]) * kx
            dy = (Pb[:, 1][:, None] - ax[:, 1][None, :]) * ky
            dd = np.hypot(dx, dy)
            km = dd.argmin(axis=1); dm = dd[np.arange(len(Pb)), km]
            sl = slice(a, a + len(Pb)); better = dm < best_d[sl]
            gi = a + np.where(better)[0]
            best_d[gi] = dm[better]; best_t[gi] = tid; best_k[gi] = km[better]
    out = []
    for i in range(N):
        if best_d[i] <= buffer:
            ax = axes[best_t[i]]; k = best_k[i]; lon, lat = ax["pts"][k]
            out.append((int(best_t[i]), round(float(lon), 4), round(float(lat), 4),
                        round(float(ax["s"][k]), 1)))
    return out

def load_axes():
    TR = {}
    for p in glob.glob(os.path.join(AXES_DIR, "*.txt")):
        m = re.match(r"(\d+)_(.+)\.txt$", os.path.basename(p))
        if not m:
            continue
        pts = [(norm360(float(ln.split()[0])), float(ln.split()[1]))
               for ln in open(p) if ln.strip() and not ln.startswith("#")]
        dens = densify(pts, AXIS_DS)
        TR[int(m.group(1))] = {"name": m.group(2).replace("_", " "),
                               "pts": np.array(dens),
                               "s": np.array(cum_dist(dens))}
    return dict(sorted(TR.items()))

def nearest_axis_point(feat_pts, axis):

    fp = np.asarray(feat_pts)
    ax = axis["pts"]
    latm = math.radians(float(ax[:, 1].mean()))
    kx = 111.320 * math.cos(latm); ky = 110.574
    dx = (fp[:, 0][:, None] - ax[:, 0][None, :]) * kx
    dy = (fp[:, 1][:, None] - ax[:, 1][None, :]) * ky
    dd = np.hypot(dx, dy)
    i, k = np.unravel_index(np.argmin(dd), dd.shape)
    return float(dd[i, k]), int(k)

def main():
    axes = load_axes()
    if not axes:
        print(f"No axes found in {AXES_DIR}/ -- run from the project folder.")
        return

    contacts = []
    missing = []
    for cls, sources in CONFIG.items():
        for path, geom in sources:
            if geom == "seamount":
                pts = read_seamounts(path, MIN_SEAMOUNT_HEIGHT_M)
                if pts is None:
                    missing.append(path); continue
                for tid, lon, lat, s_km in seamount_contacts(pts, axes, BUFFER_KM["seamount"]):
                    contacts.append((tid, cls, lon, lat, s_km))
                continue
            feats = read_features(path, geom)
            if feats is None:
                missing.append(path); continue
            feats = [densify(f, FEAT_DS) if geom != "point" and len(f) > 1 else f
                     for f in feats if f]
            for f in feats:
                for tid, axis in axes.items():
                    dist, k = nearest_axis_point(f, axis)
                    if dist <= BUFFER_KM.get(geom, 50.0):
                        lon, lat = axis["pts"][k]
                        contacts.append((tid, cls, round(float(lon), 4),
                                         round(float(lat), 4),
                                         round(float(axis["s"][k]), 1)))

    merged = []
    contacts.sort(key=lambda r: (r[0], r[1], r[4]))
    for r in contacts:
        if merged and merged[-1][0] == r[0] and merged[-1][1] == r[1] \
                and abs(merged[-1][4] - r[4]) < DEDUP_KM:
            continue
        merged.append(r)

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["trench_id", "class", "lon", "lat", "s_km"])
        w.writerows(merged)

    print(f"wrote {OUT_CSV}: {len(merged)} contacts "
          f"(from {len(contacts)} raw, {len(contacts)-len(merged)} merged)")
    by = {}
    for r in merged:
        by[r[1]] = by.get(r[1], 0) + 1
    for c in CONFIG:
        print(f"  {c:14s}: {by.get(c, 0)}")
    if missing:
        print("\n[!] These configured files were not found (edit CONFIG or add them):")
        for p in missing:
            print("    ", p)
        print("    No contacts were produced for the classes that use them.")

if __name__ == "__main__":
    main()
