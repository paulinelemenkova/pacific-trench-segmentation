#!/usr/bin/env python3

import csv, itertools, json, math, os
import numpy as np
import ruptures as rpt

HERE = os.path.dirname(os.path.abspath(__file__))
DS = 5.0
DMAX, RLMIN, WIN, FRAC, MINVALID = -2000.0, 300.0, 10, 0.8, 0.5
SMOOTH, MINSIZE, PENF, MODEL = 3, 6, 1.0, "l2"
GRID_SMOOTH, GRID_PENF, GRID_MODEL, GRID_MINSIZE = (
    (1, 3, 5),
    (0.5, 1.0, 2.0),
    ("l2", "normal"),
    (4, 6, 8),
)
STAB_TOL, STAB_MIN = 15.0, 0.6
TOLS = (15.0, 30.0, 50.0)
D_TOL = 15.0
R_PERM = 10000
SEED = 20240101
CLASSES = ("fracture_zone", "plateau", "seamount_ridge")


def load_series(tid):
    path = os.path.join(HERE, "data", "series", f"series_{tid:02d}.csv")
    with open(path) as f:
        head = f.readline()
        name = head.split(";")[0].split(" ", 3)[3]
        rows = list(csv.DictReader(f))
    S = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    S["name"] = name
    return S


def trim(S):
    v = (S["d_m"] <= DMAX) & (S["rise_land_m"] >= RLMIN)
    run = np.convolve(v, np.ones(WIN) / WIN, mode="valid")
    ok = np.where(run >= FRAC)[0]
    if not len(ok):
        return None, v
    a, b = ok[0], ok[-1] + WIN - 1
    return (a, b), v


def features(S, a, b, smooth):
    d = S["d_m"][a : b + 1]
    if smooth > 1:
        k = np.ones(smooth) / smooth
        pad = np.r_[np.full(smooth // 2, d[0]), d, np.full(smooth // 2, d[-1])]
        d = np.convolve(pad, k, mode="valid")
    gl = np.clip(S["rise_land_m"][a : b + 1], 0, None) / 20.0
    gs = np.clip(S["rise_sea_m"][a : b + 1], 0, None) / 20.0
    A = np.where(gl + gs > 0, (gl - gs) / np.maximum(gl + gs, 1e-9), 0.0)
    W = S["W_km"][a : b + 1]
    return d, W, A


def zs(x):
    sd = x.std()
    return (x - x.mean()) / (sd if sd > 0 else 1.0)


def segment(d, W, A, model=MODEL, penf=PENF, min_size=MINSIZE):
    Y = np.column_stack([zs(d), zs(W), zs(A)])
    n, p = Y.shape
    sig2 = float(np.mean(0.5 * np.var(np.diff(Y, axis=0), axis=0)))
    beta = penf * p * sig2 * math.log(n)
    if model == "normal":
        beta = penf * (p + p * (p + 1) // 2) * math.log(n)
    algo = rpt.Pelt(model=model, min_size=min_size, jump=1).fit(Y)
    bk = algo.predict(pen=beta)
    return [c for c in bk if 0 < c < n], beta, Y


def load_contacts(windows):
    C = {}
    with open(os.path.join(HERE, "data", "fabric_contacts.csv")) as f:
        for r in csv.DictReader(f):
            t = int(r["trench_id"])
            c = r["class"]
            if t not in windows:
                continue
            lo, hi = windows[t]
            s1, s2 = float(r["s1_km"]), float(r["s2_km"])
            if s2 < lo or s1 > hi:
                continue
            C.setdefault(t, {}).setdefault(c, []).append(
                (max(s1, lo) - lo, min(s2, hi) - lo, r["note"])
            )
    return C


def dist_to(b, iv):
    b = np.asarray(b, float)[:, None]
    lo = np.array([x[0] for x in iv])[None, :]
    hi = np.array([x[1] for x in iv])[None, :]
    d = np.where(b < lo, lo - b, np.where(b > hi, b - hi, 0.0))
    return d.min(axis=1)


def dist_mat(B, iv):
    lo = np.array([x[0] for x in iv])
    hi = np.array([x[1] for x in iv])
    Bx = B[..., None]
    d = np.where(Bx < lo, lo - Bx, np.where(Bx > hi, Bx - hi, 0.0))
    return d.min(axis=-1)


def rand_uniform(R, k, L, gap, rng):
    free = max(L - (k - 1) * gap, 1e-6)
    x = np.sort(rng.uniform(0, free, (R, k)), axis=1)
    return x + gap * np.arange(k)[None, :]


def coincidence(bounds, lengths, contacts, cls, tids=None, R=None, seed=SEED):
    R = R_PERM if R is None else R
    rng = np.random.default_rng(seed)
    T = [
        t
        for t in bounds
        if cls in contacts.get(t, {}) and len(bounds[t]) and (tids is None or t in tids)
    ]
    if not T:
        return None
    obs = np.concatenate([dist_to(bounds[t], contacts[t][cls]) for t in T])
    stat = lambda D: np.column_stack([D.mean(1)] + [(D <= tol).mean(1) for tol in TOLS])
    o = stat(obs[None, :])[0]
    DU, DC = [], []
    for t in T:
        b = np.asarray(bounds[t])
        L = lengths[t]
        DU.append(
            dist_mat(rand_uniform(R, len(b), L, MINSIZE * DS, rng), contacts[t][cls])
        )
        DC.append(
            dist_mat((b[None, :] + rng.uniform(0, L, (R, 1))) % L, contacts[t][cls])
        )
    nu = stat(np.concatenate(DU, axis=1))
    ns = stat(np.concatenate(DC, axis=1))

    def pv(null):
        p = [(1 + np.sum(null[:, 0] <= o[0])) / (R + 1)]
        p += [
            (1 + np.sum(null[:, j] >= o[j] - 1e-12)) / (R + 1) for j in range(1, len(o))
        ]
        return p

    return dict(
        cls=cls,
        trenches=T,
        n_b=len(obs),
        n_c=sum(len(contacts[t][cls]) for t in T),
        obs=o,
        uni_mean=nu.mean(0),
        uni_lo=np.percentile(nu, 2.5, 0),
        uni_hi=np.percentile(nu, 97.5, 0),
        p_uni=pv(nu),
        cs_mean=ns.mean(0),
        cs_lo=np.percentile(ns, 2.5, 0),
        cs_hi=np.percentile(ns, 97.5, 0),
        p_cs=pv(ns),
        null_uni=nu[:, 0],
        null_cs=ns[:, 0],
        null_uni_f=nu[:, 1:],
        null_cs_f=ns[:, 1:],
    )


CACHE_FILE = os.path.join(HERE, "results", "segmentation_cache.json")
CACHE = {}


def main():
    global CACHE
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    if os.path.exists(CACHE_FILE):
        CACHE = json.load(open(CACHE_FILE))
    summary, catalogue = [], []
    bounds, lengths, windows, series = {}, {}, {}, {}
    sens_rows = []
    for t in range(1, 21):
        S = load_series(t)
        n = len(S["s_km"])
        win, valid = trim(S)
        vfrac = valid[win[0] : win[1] + 1].mean() if win else 0.0
        row = dict(trench_id=t, trench=S["name"], axis_length_km=(n - 1) * DS)
        keep = (win[1] - win[0] + 1) / n if win else 0.0
        if win is None or vfrac < MINVALID or keep < MINVALID:
            row.update(
                status="excluded",
                reason=f"trench-floor criterion met on {valid.mean():.0%} of the axis; window keeps {keep:.0%}",
            )
            summary.append(row)
            continue
        a, b = win
        d, W, A = features(S, a, b, SMOOTH)
        cps, beta, Y = segment(d, W, A)
        s = S["s_km"][a : b + 1] - S["s_km"][a]
        bpos = [s[c] - DS / 2 for c in cps]
        bounds[t] = bpos
        lengths[t] = s[-1]
        windows[t] = (S["s_km"][a], S["s_km"][b])
        series[t] = dict(S=S, a=a, b=b, d=d, W=W, A=A, cps=cps, beta=beta)
        runs = []
        for sm, pf, md, ms in itertools.product(
            GRID_SMOOTH, GRID_PENF, GRID_MODEL, GRID_MINSIZE
        ):
            key = f"{t}|{a}|{b}|{sm}|{pf}|{md}|{ms}"
            if key not in CACHE:
                dd, WW, AA = features(S, a, b, sm)
                CACHE[key] = segment(dd, WW, AA, model=md, penf=pf, min_size=ms)[0]
            cc = CACHE[key]
            pos = np.array([s[c] - DS / 2 for c in cc])
            runs.append(pos)
            sens_rows.append(
                dict(
                    trench_id=t,
                    smooth=sm,
                    pen_factor=pf,
                    cost=md,
                    min_size=ms,
                    n_boundaries=len(cc),
                )
            )
        stab = [
            np.mean([len(r) > 0 and np.min(np.abs(r - x)) <= STAB_TOL for r in runs])
            for x in bpos
        ]
        series[t]["stab"] = stab
        seg_edges = [0] + cps + [len(d)]
        for j, c in enumerate(cps):
            l0, l1, r1 = seg_edges[j], c, seg_edges[j + 2]
            catalogue.append(
                dict(
                    trench_id=t,
                    trench=S["name"],
                    bnd_index=j + 1,
                    s_km=round(bpos[j] + S["s_km"][a], 1),
                    s_rel_km=round(bpos[j], 1),
                    lon=round(
                        0.5 * (S["pick_lon"][a + c - 1] + S["pick_lon"][a + c]), 4
                    ),
                    lat=round(
                        0.5 * (S["pick_lat"][a + c - 1] + S["pick_lat"][a + c]), 4
                    ),
                    stability=round(stab[j], 3),
                    d_left_m=round(d[l0:l1].mean()),
                    d_right_m=round(d[l1:r1].mean()),
                    W_left_km=round(W[l0:l1].mean(), 1),
                    W_right_km=round(W[l1:r1].mean(), 1),
                    A_left=round(A[l0:l1].mean(), 3),
                    A_right=round(A[l1:r1].mean(), 3),
                    dz_left=round(Y[l0:l1, 0].mean(), 3),
                    dz_right=round(Y[l1:r1, 0].mean(), 3),
                    Wz_left=round(Y[l0:l1, 1].mean(), 3),
                    Wz_right=round(Y[l1:r1, 1].mean(), 3),
                    Az_left=round(Y[l0:l1, 2].mean(), 3),
                    Az_right=round(Y[l1:r1, 2].mean(), 3),
                )
            )
        row.update(
            status="analysed",
            reason="",
            window_start_km=S["s_km"][a],
            window_end_km=S["s_km"][b],
            length_km=round(s[-1], 1),
            n_profiles=b - a + 1,
            valid_frac=round(vfrac, 3),
            depth_shallowest_m=int(S["d_m"][a : b + 1].max()),
            depth_deepest_m=int(S["d_m"][a : b + 1].min()),
            n_segments=len(cps) + 1,
            n_stable=int(np.sum(np.array(stab) >= STAB_MIN)),
        )
        summary.append(row)
        print(
            f"{t:2d} {S['name']:15s} L={s[-1]:6.0f} n_b={len(cps):3d} stable={row['n_stable']:3d}"
        )

    json.dump(CACHE, open(CACHE_FILE, "w"))
    keys = [
        "trench_id",
        "trench",
        "status",
        "reason",
        "axis_length_km",
        "window_start_km",
        "window_end_km",
        "length_km",
        "n_profiles",
        "valid_frac",
        "depth_shallowest_m",
        "depth_deepest_m",
        "n_segments",
        "n_stable",
    ]
    write("trench_summary.csv", summary, keys)
    write("boundary_catalogue.csv", catalogue, list(catalogue[0].keys()))
    write("sensitivity_grid.csv", sens_rows, list(sens_rows[0].keys()))

    contacts = load_contacts(windows)
    stable = {
        t: [x for x, st in zip(bounds[t], series[t]["stab"]) if st >= STAB_MIN]
        for t in bounds
    }
    out, bytr, loto, nulls = [], [], [], {}
    for label, B in (("all", bounds), ("stable", stable)):
        for cls in CLASSES:
            r = coincidence(B, lengths, contacts, cls)
            if r is None:
                continue
            nulls[f"{label}_{cls}"] = dict(
                obs=float(r["obs"][0]),
                uni=r["null_uni"].tolist()[:2000],
                cs=r["null_cs"].tolist()[:2000],
            )
            rec = dict(
                set=label,
                cls=cls,
                n_trenches=len(r["trenches"]),
                n_boundaries=r["n_b"],
                n_contacts=r["n_c"],
                mean_dist_km=round(r["obs"][0], 1),
                null_uni_km=round(r["uni_mean"][0], 1),
                p_uni=round(r["p_uni"][0], 4),
                null_cs_km=round(r["cs_mean"][0], 1),
                p_cs=round(r["p_cs"][0], 4),
            )
            for j, tol in enumerate(TOLS, 1):
                rec[f"frac{int(tol)}"] = round(r["obs"][j], 3)
                rec[f"frac{int(tol)}_uni"] = round(r["uni_mean"][j], 3)
                rec[f"frac{int(tol)}_uni_lo"] = round(r["uni_lo"][j], 3)
                rec[f"frac{int(tol)}_uni_hi"] = round(r["uni_hi"][j], 3)
                rec[f"p{int(tol)}_uni"] = round(r["p_uni"][j], 4)
                rec[f"frac{int(tol)}_cs"] = round(r["cs_mean"][j], 3)
                rec[f"frac{int(tol)}_cs_hi"] = round(r["cs_hi"][j], 3)
                rec[f"p{int(tol)}_cs"] = round(r["p_cs"][j], 4)
            out.append(rec)
            print(
                label,
                cls,
                {
                    k: rec[k]
                    for k in (
                        "n_trenches",
                        "n_boundaries",
                        "mean_dist_km",
                        "null_uni_km",
                        "p_uni",
                        "null_cs_km",
                        "p_cs",
                        "frac15",
                        "frac15_uni",
                        "p15_uni",
                        "p15_cs",
                    )
                },
            )
            if label == "all":
                for t in r["trenches"]:
                    rt = coincidence(B, lengths, contacts, cls, tids=[t])
                    bytr.append(
                        dict(
                            cls=cls,
                            trench_id=t,
                            trench=series[t]["S"]["name"],
                            n_boundaries=rt["n_b"],
                            n_contacts=rt["n_c"],
                            mean_dist_km=round(rt["obs"][0], 1),
                            null_uni_km=round(rt["uni_mean"][0], 1),
                            p_uni=round(rt["p_uni"][0], 4),
                            p_cs=round(rt["p_cs"][0], 4),
                            hits15=int(round(rt["obs"][1] * rt["n_b"])),
                            exp15=round(rt["uni_mean"][1] * rt["n_b"], 2),
                        )
                    )
                    if len(r["trenches"]) > 1:
                        rest = [x for x in r["trenches"] if x != t]
                        rl = coincidence(B, lengths, contacts, cls, tids=rest)
                        loto.append(
                            dict(
                                cls=cls,
                                left_out=series[t]["S"]["name"],
                                n_boundaries=rl["n_b"],
                                mean_dist_km=round(rl["obs"][0], 1),
                                null_uni_km=round(rl["uni_mean"][0], 1),
                                p_uni=round(rl["p_uni"][0], 4),
                                p_cs=round(rl["p_cs"][0], 4),
                                frac15=round(rl["obs"][1], 3),
                                p15_uni=round(rl["p_uni"][1], 4),
                            )
                        )
    write("coincidence_by_class.csv", out, list(out[0].keys()))
    write("coincidence_by_trench.csv", bytr, list(bytr[0].keys()))
    write("coincidence_loto.csv", loto, list(loto[0].keys()))
    json.dump(
        nulls, open(os.path.join(HERE, "results", "null_distributions.json"), "w")
    )
    rows = []
    for t, cc in contacts.items():
        for cls, iv in cc.items():
            for s1, s2, note in iv:
                rows.append(
                    dict(
                        trench_id=t,
                        cls=cls,
                        s1_rel_km=s1,
                        s2_rel_km=s2,
                        s1_km=s1 + windows[t][0],
                        s2_km=s2 + windows[t][0],
                        note=note,
                    )
                )
    write("contacts_in_windows.csv", rows, list(rows[0].keys()))


def write(name, rows, keys):
    with open(os.path.join(HERE, "results", name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


if __name__ == "__main__":
    main()
