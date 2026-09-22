#!/usr/bin/env python3
import glob
import math
import os
import re
import csv
import statistics

import numpy as np
import pandas as pd
import pygmt

GRID       = "@earth_relief_01m"
DS         = 5.0
L          = 60.0
DU         = 3.0
SMOOTH     = 3
LMIN       = 6
D_TOL      = 15.0
R_PERM     = 10000
SEED       = 20240101
AXES_DIR   = "axes_full"
FABRIC_CSV = "fabric_intersections.csv"
R_EARTH    = 6371.0
CLASSES    = ["fracture_zone", "plateau", "ridge"]

def hav(a, b):
    lo1, la1 = map(math.radians, a); lo2, la2 = map(math.radians, b)
    h = math.sin((la2 - la1) / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin((lo2 - lo1) / 2) ** 2
    return 2 * R_EARTH * math.asin(math.sqrt(h))

def bearing(a, b):
    lo1, la1 = map(math.radians, a); lo2, la2 = map(math.radians, b)
    y = math.sin(lo2 - lo1) * math.cos(la2)
    x = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(lo2 - lo1)
    return math.atan2(y, x)

def dest(lon, lat, brg, dist):
    la1, lo1, dR = math.radians(lat), math.radians(lon), dist / R_EARTH
    la2 = math.asin(math.sin(la1) * math.cos(dR) + math.cos(la1) * math.sin(dR) * math.cos(brg))
    lo2 = lo1 + math.atan2(math.sin(brg) * math.sin(dR) * math.cos(la1),
                           math.cos(dR) - math.sin(la1) * math.sin(la2))
    return (math.degrees(lo2) % 360, math.degrees(la2))

def densify(ax):
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
    return dense

def width_asymmetry(profile, uf, zf):

    thr = zf + 2000.0
    kf = min(range(len(profile)), key=lambda k: abs(profile[k][0] - uf))
    kL = kf
    while kL > 0 and profile[kL][1] < thr:
        kL -= 1
    kR = kf
    while kR < len(profile) - 1 and profile[kR][1] < thr:
        kR += 1
    uL, zL = profile[kL]; uR, zR = profile[kR]
    W = uR - uL
    gL = (zL - zf) / max(uf - uL, 1e-6)
    gR = (zR - zf) / max(uR - uf, 1e-6)
    A = (gL - gR) / max(gL + gR, 1e-6)
    return W, A

def sample_trench(dense):

    n = len(dense)
    brg = [bearing(dense[max(i - 1, 0)], dense[min(i + 1, n - 1)]) for i in range(n)]
    offs = [-L + k * DU for k in range(int(2 * L // DU) + 1)]
    lons, lats, idx, us = [], [], [], []
    for i, (lon, lat) in enumerate(dense):
        for u in offs:
            b = brg[i] + math.pi / 2 if u >= 0 else brg[i] - math.pi / 2
            plon, plat = dest(lon, lat, b, abs(u))
            lons.append(plon); lats.append(plat); idx.append(i); us.append(u)
    z = pygmt.grdtrack(points=pd.DataFrame({"x": lons, "y": lats}),
                       grid=GRID, newcolname="z")["z"].to_list()
    prof = {}
    for i, u, zz in zip(idx, us, z):
        if zz == zz:
            prof.setdefault(i, []).append((u, zz))
    s, lonc, latc, d, W, A = [], [], [], [], [], []
    for i in range(n):
        pr = prof.get(i)
        if not pr or len(pr) < 5:
            continue
        pr.sort()
        cand = [(zz, k) for k, (u, zz) in enumerate(pr) if abs(u) <= 50]
        if not cand:
            continue
        zf, kf = min(cand); uf = pr[kf][0]
        w, a = width_asymmetry(pr, uf, zf)
        s.append(i * DS); lonc.append(dense[i][0]); latc.append(dense[i][1])
        d.append(zf); W.append(w); A.append(a)
    return s, lonc, latc, d, W, A

def zscore(a):
    mu = sum(a) / len(a); sd = (sum((x - mu) ** 2 for x in a) / len(a)) ** 0.5 or 1.0
    return [(x - mu) / sd for x in a]

def segment(d, W, A):

    Y = np.array([zscore(d), zscore(W), zscore(A)]).T
    n, p = Y.shape
    sig2 = sum(0.5 * statistics.pvariance(Y[1:, j] - Y[:-1, j]) for j in range(p)) / p
    beta = p * sig2 * math.log(n)
    try:
        import ruptures as rpt
        algo = rpt.Pelt(model="rbf", min_size=LMIN).fit(Y)
        bkps = algo.predict(pen=beta)
        return [b for b in bkps if 0 < b < n]
    except Exception:
        return _optimal_partition(Y, beta)

def _optimal_partition(Y, beta):
    n, p = Y.shape
    cs1 = np.vstack([np.zeros(p), np.cumsum(Y, axis=0)])
    cs2 = np.vstack([np.zeros(p), np.cumsum(Y ** 2, axis=0)])

    def sc(a, b):
        nn = b - a
        return float(np.sum(cs2[b] - cs2[a] - (cs1[b] - cs1[a]) ** 2 / nn))

    F = [None] * (n + 1); F[0] = -beta; prev = [0] * (n + 1)
    for t in range(1, n + 1):
        hi = t - LMIN
        cands = range(0, hi + 1) if hi >= 0 else [0]
        best = None; bs = 0
        for a in cands:
            c = F[a] + sc(a, t) + beta
            if best is None or c < best:
                best = c; bs = a
        F[t] = best; prev[t] = bs
    cps = []; t = n
    while t > 0:
        a = prev[t]
        if a > 0:
            cps.append(a)
        t = a
    return sorted(cps)

def load_axes():
    TR = {}
    for p in glob.glob(os.path.join(AXES_DIR, "*.txt")):
        m = re.match(r"(\d+)_(.+)\.txt$", os.path.basename(p))
        if not m:
            continue
        TR[int(m.group(1))] = (m.group(2).replace("_", " "),
            [tuple(map(float, ln.split()[:2])) for ln in open(p)
             if ln.strip() and not ln.startswith("#")])
    return dict(sorted(TR.items()))

def nearest(a, arr):
    return min((abs(a - x) for x in arr), default=float("inf"))

def coincidence(per_trench_bounds, per_trench_contacts, lengths):

    rng = np.random.default_rng(SEED)
    rows = []
    for cls in CLASSES:
        dobs, hits, nb, nc = [], 0, 0, 0
        for tid, bnds in per_trench_bounds.items():
            con = per_trench_contacts.get(tid, {}).get(cls, [])
            nc += len(con)
            if not con:
                continue
            for b in bnds:
                dist = nearest(b, con); dobs.append(dist)
                hits += (dist <= D_TOL); nb += 1
        if not dobs:
            rows.append((cls, 0, nc, float("nan"), float("nan"), float("nan"), float("nan")))
            continue
        D_obs = float(np.mean(dobs)); frac = hits / nb

        null = np.empty(R_PERM)
        min_gap = LMIN * DS
        for r in range(R_PERM):
            dd = []
            for tid, bnds in per_trench_bounds.items():
                con = per_trench_contacts.get(tid, {}).get(cls, [])
                if not con or not bnds:
                    continue
                rb = _random_boundaries(len(bnds), lengths[tid], min_gap, rng)
                dd += [nearest(b, con) for b in rb]
            null[r] = np.mean(dd) if dd else np.nan
        null = null[~np.isnan(null)]
        p = (1 + int(np.sum(null <= D_obs))) / (len(null) + 1)
        rows.append((cls, nb, nc, D_obs, float(np.mean(null)), frac, p))
    return rows

def _random_boundaries(k, length, min_gap, rng):
    for _ in range(200):
        xs = np.sort(rng.uniform(0, length, k))
        if k < 2 or np.all(np.diff(xs) >= min_gap):
            return xs.tolist()
    return np.sort(rng.uniform(0, length, k)).tolist()

def main():
    TR = load_axes()
    summary, catalogue = [], []
    per_bounds, lengths, node_lookup = {}, {}, {}

    for tid, (name, ax) in TR.items():
        dense = densify(ax)
        s, lon, lat, d, W, A = sample_trench(dense)
        if len(d) < 2 * LMIN:
            print(f"  ! trench {tid} {name}: too few nodes, skipped"); continue
        cps = segment(d, W, A)
        bpos = [s[c] for c in cps]
        per_bounds[tid] = bpos
        lengths[tid] = s[-1]
        node_lookup[tid] = (s, lon, lat)
        summary.append({"trench_id": tid, "trench": name, "length_km": round(s[-1], 1),
                        "n_profiles": len(s), "depth_min_m": int(max(d)),
                        "depth_max_m": int(min(d)), "n_segments": len(cps) + 1})
        for j, c in enumerate(cps, 1):
            catalogue.append({"trench_id": tid, "trench": name, "bnd_index": j,
                              "s_km": round(s[c], 1), "lon": round(lon[c], 4),
                              "lat": round(lat[c], 4)})
        print(f"  trench {tid:2d} {name:16s}: {len(s)} profiles, {len(cps)} boundaries")

    _write("trench_summary.csv", summary,
           ["trench_id", "trench", "length_km", "n_profiles",
            "depth_min_m", "depth_max_m", "n_segments"])
    _write("boundary_catalogue.csv", catalogue,
           ["trench_id", "trench", "bnd_index", "s_km", "lon", "lat"])
    print("wrote trench_summary.csv and boundary_catalogue.csv")

    if not os.path.exists(FABRIC_CSV):
        print(f"\n[skip] {FABRIC_CSV} not found -- see RESULTS_HOWTO.md, Part B.")
        return
    contacts = _load_contacts(node_lookup)
    rows = coincidence(per_bounds, contacts, lengths)
    _write("coincidence_by_class.csv",
           [dict(zip(["class", "n_boundaries", "n_contacts", "mean_dist_km",
                      "null_mean_km", "fraction", "p_value"],
                     (c, nb, nc, round(md, 2), round(nm, 2), round(fr, 3), round(pv, 4))))
            for (c, nb, nc, md, nm, fr, pv) in rows],
           ["class", "n_boundaries", "n_contacts", "mean_dist_km",
            "null_mean_km", "fraction", "p_value"])
    print("wrote coincidence_by_class.csv")

def _load_contacts(node_lookup):

    contacts = {}
    with open(FABRIC_CSV) as f:
        for row in csv.DictReader(f):
            tid = int(row["trench_id"]); cls = row["class"].strip()
            if cls not in CLASSES or tid not in node_lookup:
                continue
            if row.get("s_km") not in (None, ""):
                s_km = float(row["s_km"])
            else:
                s, lon, lat = node_lookup[tid]
                lo, la = float(row["lon"]) % 360, float(row["lat"])
                k = min(range(len(s)),
                        key=lambda i: (lon[i] - lo) ** 2 + (lat[i] - la) ** 2)
                s_km = s[k]
            contacts.setdefault(tid, {}).setdefault(cls, []).append(s_km)
    return contacts

def _write(path, rows, header):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header); w.writeheader()
        for r in rows:
            w.writerow(r)

if __name__ == "__main__":
    main()
