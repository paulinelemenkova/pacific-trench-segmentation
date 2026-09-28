#!/usr/bin/env python3

import base64, csv, hashlib, math, subprocess, sys, zlib
import numpy as np
from shapely.geometry import shape, Point, Polygon
from shapely.ops import nearest_points

R = 6371.0
DS = 5.0
D_FZ, X_FZ, D_PL, D_SM, H_SM = 30.0, 200.0, 50.0, 60.0, 1000.0
SIGN = {
    1: 1,
    2: 1,
    3: -1,
    4: 1,
    5: 1,
    6: 1,
    7: -1,
    8: 1,
    9: 1,
    10: -1,
    11: 1,
    12: 1,
    13: 1,
    14: 1,
    15: -1,
    16: -1,
    17: 1,
    18: -1,
    19: -1,
    20: -1,
}


def read_axes(path):
    T = []
    cur = None
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        if ln.startswith(">"):
            cur = [ln.split('"')[1], []]
            T.append(cur)
            continue
        x, y = map(float, ln.split()[:2])
        cur[1].append((x % 360, y))
    return [(n, np.array(p)) for n, p in T]


def havv(lo1, la1, lo2, la2):
    lo1, la1, lo2, la2 = map(np.radians, (lo1, la1, lo2, la2))
    h = (
        np.sin((la2 - la1) / 2) ** 2
        + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    )
    return 2 * R * np.arcsin(np.sqrt(np.clip(h, 0, 1)))


def bearing(a, b):
    lo1, la1, lo2, la2 = map(math.radians, (*a, *b))
    return math.atan2(
        math.sin(lo2 - lo1) * math.cos(la2),
        math.cos(la1) * math.sin(la2)
        - math.sin(la1) * math.cos(la2) * math.cos(lo2 - lo1),
    )


def densify_axis(ax):
    d = havv(ax[:-1, 0], ax[:-1, 1], ax[1:, 0], ax[1:, 1])
    cum = np.r_[0, np.cumsum(d)]
    out = []
    for s in np.arange(0, cum[-1] + 1e-9, DS):
        j = min(np.searchsorted(cum, s, side="right") - 1, len(ax) - 2)
        t = (s - cum[j]) / max(cum[j + 1] - cum[j], 1e-9)
        out.append(ax[j] + t * (ax[j + 1] - ax[j]))
    return np.array(out)


def densify_line(p, step=5.0):
    out = [p[0]]
    for a, b in zip(p[:-1], p[1:]):
        k = max(int(havv(a[0], a[1], b[0], b[1]) // step), 1)
        d = np.array([(b[0] - a[0] + 180) % 360 - 180, b[1] - a[1]])
        for m in range(1, k + 1):
            q = a + d * m / k
            out.append(np.array([q[0] % 360, q[1]]))
    return np.array(out)


class Axis:
    def __init__(self, nodes, sign):
        self.p = nodes
        n = len(nodes)
        self.s = np.arange(n) * DS
        self.brg = np.array(
            [bearing(nodes[max(i - 2, 0)], nodes[min(i + 2, n - 1)]) for i in range(n)]
        )
        self.sign = sign

    def project(self, lon, lat):
        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)
        D = havv(lon[:, None], lat[:, None], self.p[None, :, 0], self.p[None, :, 1])
        k = D.argmin(1)
        dmin = D[np.arange(len(k)), k]
        la0 = np.radians(self.p[k, 1])
        dx = ((lon - self.p[k, 0] + 180) % 360 - 180) * np.cos(la0) * 111.32
        dy = (lat - self.p[k, 1]) * 110.57
        b = self.brg[k]
        along = dx * np.sin(b) + dy * np.cos(b)
        cross = (dx * np.cos(b) - dy * np.sin(b)) * self.sign
        return k, dmin, cross, along


def read_gmt_lines(path):
    feats, cur = [], []
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        if ln.startswith(">"):
            if cur:
                feats.append(np.array(cur))
                cur = []
            continue
        x, y = map(float, ln.split()[:2])
        cur.append((x % 360, y))
    if cur:
        feats.append(np.array(cur))
    return feats


def main():
    axes = read_axes(sys.argv[1])
    subprocess.run(
        ["ogr2ogr", "-f", "GMT", "fz.gmt", "SHP/GSFML_SF_FZ_KM.shp"], check=True
    )
    fzs = [densify_line(f) for f in read_gmt_lines("fz.gmt") if len(f) > 1]
    import json

    subprocess.run(
        [
            "ogr2ogr",
            "-f",
            "GeoJSON",
            "lips.json",
            "IgneousProvinces/Whittaker_2015/SHP/Whittaker_etal_2015_LIPs.shp",
        ],
        check=True,
    )
    lips = []
    if True:
        for rec in json.load(open("lips.json"))["features"]:
            g = shape(rec["geometry"])
            nm = rec["properties"]["NAME"] or "LIP"
            polys = [g] if isinstance(g, Polygon) else list(g.geoms)
            for pg in polys:
                xy = np.array(pg.exterior.coords)
                xy[:, 0] %= 360
                if (
                    np.ptp(xy[:, 0]) > 180
                ):
                    continue
                lips.append((nm, Polygon(xy)))
    sm = []
    for ln in open("kw.txt"):
        if ln.startswith("#") or ln.startswith(">") or not ln.strip():
            continue
        p = ln.split()
        if float(p[5]) >= H_SM:
            sm.append((float(p[0]) % 360, float(p[1]), float(p[5])))
    sm = np.array(sm)
    rows = []
    geo_fz = []
    geo_pl = {}
    for tid, (name, ax) in enumerate(axes, 1):
        A = Axis(densify_axis(ax), SIGN[tid])
        P = A.p
        w, e, s_, n_ = (
            P[:, 0].min() - 4,
            P[:, 0].max() + 4,
            P[:, 1].min() - 4,
            P[:, 1].max() + 4,
        )
        for fi, f in enumerate(fzs):
            m = (f[:, 0] > w) & (f[:, 0] < e) & (f[:, 1] > s_) & (f[:, 1] < n_)
            if m.sum() < 2:
                continue
            g = f[m]
            k, dmin, cr, al = A.project(g[:, 0], g[:, 1])
            j = dmin.argmin()
            if dmin[j] <= D_FZ and 0 < k[j] < len(P) - 1:
                rows.append(
                    (
                        tid,
                        "fracture_zone",
                        A.s[k[j]],
                        A.s[k[j]],
                        P[k[j], 0],
                        P[k[j], 1],
                        "direct",
                    )
                )
                geo_fz.append(g)
                continue
            for end, prev in (
                (g[0], g[min(10, len(g) - 1)]),
                (g[-1], g[max(-11, -len(g))]),
            ):
                ke, de, ce, ae = A.project(end[0], end[1])
                if de[0] > X_FZ or ce[0] < 0:
                    continue
                v = end - prev
                nv = math.hypot(v[0] * math.cos(math.radians(end[1])), v[1])
                if nv == 0:
                    continue
                hit = None
                for step in np.arange(5.0, X_FZ + 5.0, 5.0):
                    q = end + v * (step / (nv * 111.0))
                    kq, dq, cq, aq = A.project(q[0], q[1])
                    if cq[0] <= 0 and 0 < kq[0] < len(P) - 1:
                        hit = kq[0]
                        break
                if hit is not None:
                    rows.append(
                        (
                            tid,
                            "fracture_zone",
                            A.s[hit],
                            A.s[hit],
                            P[hit, 0],
                            P[hit, 1],
                            "extrap%.0f" % step,
                        )
                    )
                    geo_fz.append(g)
                    break
        for li, (nm, pg) in enumerate(lips):
            bx = pg.bounds
            if bx[2] < w or bx[0] > e or bx[3] < s_ or bx[1] > n_:
                continue
            hitk = []
            for i, (lo, la) in enumerate(P):
                pt = Point(lo, la)
                if pg.contains(pt):
                    hitk.append(i)
                    continue
                q = nearest_points(pg, pt)[0]
                if havv(lo, la, q.x, q.y) <= D_PL:
                    hitk.append(i)
            if not hitk:
                continue
            runs = []
            st = hitk[0]
            pr = hitk[0]
            for i in hitk[1:]:
                if i - pr > 2:
                    runs.append((st, pr))
                    st = i
                pr = i
            runs.append((st, pr))
            for a, b in runs:
                mid = (a + b) // 2
                rows.append(
                    (
                        tid,
                        "plateau",
                        A.s[a],
                        A.s[b],
                        P[mid, 0],
                        P[mid, 1],
                        nm.replace(",", " "),
                    )
                )
            geo_pl[li] = (nm, np.array(pg.exterior.coords))
        m = (sm[:, 0] > w) & (sm[:, 0] < e) & (sm[:, 1] > s_) & (sm[:, 1] < n_)
        if m.any():
            g = sm[m]
            k, dmin, cr, al = A.project(g[:, 0], g[:, 1])
            ok = (
                (cr >= -10)
                & (cr <= D_SM)
                & (np.abs(al) <= 10)
                & (k > 0)
                & (k < len(P) - 1)
            )
            for i in np.where(ok)[0]:
                rows.append(
                    (
                        tid,
                        "seamount_ridge",
                        A.s[k[i]],
                        A.s[k[i]],
                        g[i, 0],
                        g[i, 1],
                        "h%.0f" % g[i, 2],
                    )
                )
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    out = []
    for r in rows:
        if out and out[-1][0] == r[0] and out[-1][1] == r[1]:
            if r[1] != "plateau" and abs(out[-1][2] - r[2]) < 25:
                continue
            if r[1] == "plateau" and out[-1][6] == r[6] and r[2] <= out[-1][3] + 10:
                o = out[-1]
                out[-1] = (o[0], o[1], o[2], max(o[3], r[3]), o[4], o[5], o[6])
                continue
        out.append(r)
    with open("contacts.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["trench_id", "class", "s1_km", "s2_km", "lon", "lat", "note"])
        for r in out:
            wr.writerow(
                [
                    r[0],
                    r[1],
                    "%.0f" % r[2],
                    "%.0f" % r[3],
                    "%.3f" % r[4],
                    "%.3f" % r[5],
                    r[6],
                ]
            )
    print(
        "contacts",
        len(out),
        {
            c: sum(1 for r in out if r[1] == c)
            for c in ("fracture_zone", "plateau", "seamount_ridge")
        },
    )
    txt = open("contacts.csv").read()
    raw = zlib.compress(txt.encode(), 9)
    s = base64.b64encode(raw).decode()
    print(f"B 101 0 {hashlib.md5(raw).hexdigest()[:8]} {len(s)}")
    for k in range(0, len(s), 100):
        print(s[k : k + 100])
    lines = []
    for g in geo_fz:
        lines.append(("F", g[::4]))
    for li, (nm, xy) in geo_pl.items():
        xy = densify_line(np.c_[xy[:, 0] % 360, xy[:, 1]], 20.0)
        lines.append(("P", xy[::1]))
    buf = []
    for t, xy in lines:
        buf.append(f"> {t}")
        buf += ["%.2f %.2f" % (x, y) for x, y in xy]
    txt = "\n".join(buf) + "\n"
    raw = zlib.compress(txt.encode(), 9)
    s = base64.b64encode(raw).decode()
    print(f"B 102 0 {hashlib.md5(raw).hexdigest()[:8]} {len(s)}")
    for k in range(0, len(s), 100):
        print(s[k : k + 100])


if __name__ == "__main__":
    main()
