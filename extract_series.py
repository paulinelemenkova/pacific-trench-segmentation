#!/usr/bin/env python3

import base64, hashlib, math, subprocess, sys, zlib
import numpy as np

GRID = "@earth_gebco_30s"
DS = 5.0
DU = 2.5
UMAX = 120.0
R = 6371.0
UPICK = 60.0
LAMBDA = 60.0
HW = 1000.0
WCAP = 60.0


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


def hav(a, b):
    lo1, la1, lo2, la2 = map(math.radians, (*a, *b))
    h = (
        math.sin((la2 - la1) / 2) ** 2
        + math.cos(la1) * math.cos(la2) * math.sin((lo2 - lo1) / 2) ** 2
    )
    return 2 * R * math.asin(math.sqrt(h))


def bearing(a, b):
    lo1, la1, lo2, la2 = map(math.radians, (*a, *b))
    return math.atan2(
        math.sin(lo2 - lo1) * math.cos(la2),
        math.cos(la1) * math.sin(la2)
        - math.sin(la1) * math.cos(la2) * math.cos(lo2 - lo1),
    )


def dest(lon, lat, brg, d):
    la1, lo1, dr = math.radians(lat), math.radians(lon), d / R
    la2 = math.asin(
        math.sin(la1) * math.cos(dr) + math.cos(la1) * math.sin(dr) * math.cos(brg)
    )
    lo2 = lo1 + math.atan2(
        math.sin(brg) * math.sin(dr) * math.cos(la1),
        math.cos(dr) - math.sin(la1) * math.sin(la2),
    )
    return math.degrees(lo2) % 360, math.degrees(la2)


def densify(ax):
    cum = [0.0]
    for i in range(1, len(ax)):
        cum.append(cum[-1] + hav(ax[i - 1], ax[i]))
    cum = np.array(cum)
    out = []
    for s in np.arange(0, cum[-1] + 1e-9, DS):
        j = min(np.searchsorted(cum, s, side="right") - 1, len(ax) - 2)
        t = (s - cum[j]) / max(cum[j + 1] - cum[j], 1e-9)
        out.append(ax[j] + t * (ax[j + 1] - ax[j]))
    return np.array(out)


def profiles(d, u):
    n = len(d)
    brg = np.array([bearing(d[max(i - 2, 0)], d[min(i + 2, n - 1)]) for i in range(n)])
    pts = [
        dest(
            d[i, 0],
            d[i, 1],
            brg[i] + (math.pi / 2 if uu >= 0 else -math.pi / 2),
            abs(uu),
        )
        for i in range(n)
        for uu in u
    ]
    P = np.array(pts)
    np.savetxt("pts.txt", P, fmt="%.5f")
    reg = "-R%.2f/%.2f/%.2f/%.2f" % (
        P[:, 0].min() - 0.2,
        P[:, 0].max() + 0.2,
        P[:, 1].min() - 0.2,
        P[:, 1].max() + 0.2,
    )
    subprocess.run(
        ["gmt", "grdcut", GRID, reg, "-Gcut.nc"], check=True, capture_output=True
    )
    out = subprocess.run(
        ["gmt", "grdtrack", "pts.txt", "-Gcut.nc", "-N"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    z = np.array([float(l.split()[2]) for l in out.splitlines() if l.strip()])
    return z.reshape(n, len(u))


def viterbi(Zw, step):
    n, m = Zw.shape
    Zf = np.where(np.isnan(Zw), 0.0, Zw)
    idx = np.arange(m)
    trans = LAMBDA * step * np.abs(idx[:, None] - idx[None, :])
    cost = Zf[0].copy()
    back = np.zeros((n, m), int)
    for i in range(1, n):
        tot = cost[:, None] + trans
        back[i] = np.argmin(tot, axis=0)
        cost = tot[back[i], idx] + Zf[i]
    path = np.empty(n, int)
    path[-1] = int(np.argmin(cost))
    for i in range(n - 1, 0, -1):
        path[i - 1] = back[i, path[i]]
    return path


def metrics(Z, u, sign):
    Z = Z[:, ::-1] if sign < 0 else Z
    w = np.abs(u) <= UPICK
    uw = u[w]
    path = viterbi(Z[:, w], DU)
    up = uw[path]
    n = len(Z)
    d = np.array([Z[i, np.searchsorted(u, up[i])] for i in range(n)])

    def at(i, x):
        return np.interp(x, u, Z[i])

    rl = np.array([at(i, up[i] - 20) - d[i] for i in range(n)])
    rs = np.array([at(i, up[i] + 20) - d[i] for i in range(n)])
    W = np.empty(n)
    for i in range(n):
        k0 = np.searchsorted(u, up[i])
        thr = d[i] + HW
        kL = k0
        while kL > 0 and Z[i, kL] < thr and u[k0] - u[kL] < WCAP:
            kL -= 1
        kR = k0
        while kR < len(u) - 1 and Z[i, kR] < thr and u[kR] - u[k0] < WCAP:
            kR += 1
        W[i] = u[kR] - u[kL]
    Ro = np.array(
        [
            np.nanmedian(Z[i, (u >= up[i] + 30) & (u <= up[i] + 60)]) - d[i]
            for i in range(n)
        ]
    )
    return up, d, rl, rs, W, Ro


def main():
    args = sys.argv[1:]
    b64 = "--b64" in args
    ids = [int(a) for a in args[1:] if a != "--b64"]
    axes = read_axes(args[0])
    u = np.arange(-UMAX, UMAX + 1e-9, DU)
    for tid in ids:
        name, ax = axes[tid - 1]
        dn = densify(ax)
        Z = profiles(dn, u)
        far = np.nanmean(Z[:, u >= 60]) - np.nanmean(Z[:, u <= -60])
        sign = 1 if far < 0 else -1
        up, d, rl, rs, W, Ro = metrics(Z, u, sign)
        n = len(d)
        with open(f"series_{tid:02d}.csv", "w") as f:
            f.write("s_km,u_pick_km,d_m,rise_land_m,rise_sea_m,W_km,relief_outer_m\n")
            for i in range(n):
                f.write(
                    f"{i * DS:.1f},{up[i]:.1f},{d[i]:.0f},{rl[i]:.0f},{rs[i]:.0f},{W[i]:.1f},{Ro[i]:.0f}\n"
                )
        print(
            f"T {tid} {name.replace(' ', '_')} n={n} sign={sign} dmin={d.min():.0f} dmax={d.max():.0f}",
            flush=True,
        )
        if b64:
            q = np.vstack(
                [
                    np.round(up / DU),
                    np.round(d / 5),
                    np.round(rl / 10),
                    np.round(rs / 10),
                    np.round(W / DU),
                    np.round(Ro / 10),
                ]
            ).astype(np.int32)
            q = np.diff(q, axis=1, prepend=0).astype(np.int16)
            raw = zlib.compress(q.tobytes(), 9)
            s = base64.b64encode(raw).decode()
            print(f"B {tid} {n} {hashlib.md5(raw).hexdigest()[:8]} {len(s)}")
            for k in range(0, len(s), 100):
                print(s[k : k + 100])


if __name__ == "__main__":
    main()
