#!/usr/bin/env python3
import csv
import numpy as np
import pygmt

DTOL = 15.0
MIN_GAP = 30.0
R = 10000
SEED = 20240101
CLASSES = ["fracture_zone", "plateau", "ridge"]
LABEL = {"fracture_zone": "Fracture zones", "plateau": "Oceanic plateaus",
         "ridge": "Seamounts / ridges"}
HISTCOL = {"fracture_zone": "steelblue", "plateau": "seagreen", "ridge": "mediumpurple"}
BARCOL = {"fracture_zone": "firebrick", "plateau": "forestgreen", "ridge": "magenta"}

def load():
    bnd = {}
    for r in csv.DictReader(open("boundary_catalogue.csv")):
        bnd.setdefault(int(r["trench_id"]), []).append(float(r["s_km"]))
    con = {}
    for r in csv.DictReader(open("fabric_intersections.csv")):
        con.setdefault(int(r["trench_id"]), {}).setdefault(r["class"], []).append(float(r["s_km"]))
    length = {int(r["trench_id"]): float(r["length_km"])
              for r in csv.DictReader(open("trench_summary.csv"))}
    return bnd, con, length

def run_test(bnd, con, length):
    rng = np.random.default_rng(SEED)

    def randb(k, L):
        for _ in range(200):
            xs = np.sort(rng.uniform(0, L, k))
            if k < 2 or np.all(np.diff(xs) >= MIN_GAP):
                return xs
        return np.sort(rng.uniform(0, L, k))

    out = {}
    for c in CLASSES:
        per, dobs = [], []
        for tid, bl in bnd.items():
            cc = con.get(tid, {}).get(c, [])
            if not cc:
                continue
            per.append((tid, bl, np.array(cc)))
            for b in bl:
                dobs.append(min(abs(b - f) for f in cc))
        dobs = np.array(dobs)
        nulls = np.empty(R)
        nullfrac = np.empty(R)
        for i in range(R):
            dd = []
            for tid, bl, cc in per:
                for b in randb(len(bl), length[tid]):
                    dd.append(np.min(np.abs(b - cc)))
            dd = np.array(dd)
            nulls[i] = dd.mean()
            nullfrac[i] = (dd <= DTOL).mean()
        out[c] = dict(nb=dobs.size, obs=dobs.mean(), nulls=nulls,
                      p=(1 + int(np.sum(nulls <= dobs.mean()))) / (R + 1),
                      fobs=float((dobs <= DTOL).mean()),
                      fnull=float(nullfrac.mean()),
                      flo=float(np.percentile(nullfrac, 2.5)),
                      fhi=float(np.percentile(nullfrac, 97.5)))
    return out

def fig_coincidence(res):
    fig = pygmt.Figure()
    pygmt.config(FONT_HEADING="13p,Helvetica-Bold", FONT_TITLE="11p,Helvetica-Bold",
                 FONT_ANNOT_PRIMARY="8p", FONT_LABEL="9p", MAP_FRAME_TYPE="plain",
                 MAP_GRID_PEN_PRIMARY="0.2p,gray88")
    tags = ["(a)", "(b)", "(c)"]
    with fig.subplot(nrows=1, ncols=3, figsize=("16.2c", "5.2c"), margins=["0.85c", "0.5c"],
                     title="Permutation null distributions of mean boundary-to-fabric distance"):
        for j, c in enumerate(CLASSES):
            d = res[c]
            nulls, obs = d["nulls"], d["obs"]
            lo = min(nulls.min(), obs); hi = max(nulls.max(), obs)
            pad = 0.06 * (hi - lo); lo -= pad; hi += pad
            cnt, edges = np.histogram(nulls, bins=40, range=(lo, hi))
            yt = int(cnt.max() * 1.2)
            with fig.set_panel(j):
                fig.basemap(region=[lo, hi, 0, yt], projection="X5.4c/5.2c",
                            frame=["WSne+t%s %s" % (tags[j], LABEL[c]),
                                   "xaf+lMean distance (km)",
                                   "yaf" + ("+lNull frequency" if j == 0 else "")])
                for k in range(len(cnt)):
                    if cnt[k] > 0:
                        fig.plot(x=[edges[k], edges[k], edges[k+1], edges[k+1]],
                                 y=[0, cnt[k], cnt[k], 0], fill=HISTCOL[c],
                                 pen="0.2p,gray30", close=True)
                fig.plot(x=[obs, obs], y=[0, yt], pen="1.4p,red")
                fig.text(x=obs, y=yt, text="observed", font="8p,Helvetica-Bold,red",
                         justify="RT" if j != 1 else "LT",
                         offset=("-0.1c/-0.15c" if j != 1 else "0.1c/-0.15c"), no_clip=True)
                px = lo if j != 1 else hi
                fig.text(x=px, y=yt, text="p = %.4f" % d["p"], font="8p,Helvetica,black",
                         justify="LT" if j != 1 else "RT",
                         offset=("0.15c/-0.15c" if j != 1 else "-0.15c/-0.15c"), no_clip=True)
    fig.savefig("fig_coincidence.png", dpi=300)
    fig.savefig("fig_coincidence.pdf")

def fig_classfrac(res):
    fig = pygmt.Figure()
    pygmt.config(FONT_TITLE="12p,Helvetica-Bold", FONT_ANNOT_PRIMARY="9p", FONT_LABEL="10p",
                 MAP_FRAME_TYPE="plain", MAP_GRID_PEN_PRIMARY="0.2p,gray88")
    ymax = max(res[c]["fhi"] for c in CLASSES) * 1.3
    fig.basemap(region=[0.4, 3.6, 0, ymax], projection="X11c/7.2c",
                frame=["WSne+tCoincident-boundary fraction: observed vs permutation null",
                       "yaf+lFraction of boundaries within @~\\243@~15 km"])
    for i, c in enumerate(CLASSES, 1):
        d = res[c]
        fig.plot(x=[i-0.32, i-0.32, i+0.32, i+0.32], y=[0, d["fobs"], d["fobs"], 0],
                 fill=BARCOL[c] + "@25", pen="0.5p,black", close=True)
        fig.plot(x=[i, i], y=[d["flo"], d["fhi"]], pen="1.4p,gray20")
        fig.plot(x=[i-0.12, i+0.12], y=[d["flo"], d["flo"]], pen="1.4p,gray20")
        fig.plot(x=[i-0.12, i+0.12], y=[d["fhi"], d["fhi"]], pen="1.4p,gray20")
        fig.plot(x=[i], y=[d["fnull"]], style="d0.28c", fill="white", pen="1.2p,gray20")
        pv = "p<0.001" if d["p"] < 0.001 else "p=%.2f" % d["p"]
        fig.text(x=i, y=d["fhi"] + 0.04 * ymax, text=pv,
                 font="8p,Helvetica-Oblique,gray25", justify="CB", no_clip=True)
    for i, c in enumerate(CLASSES, 1):
        fig.text(x=i, y=-0.05 * ymax, text=LABEL[c].replace(" / ", " /@@"),
                 font="9p,Helvetica-Bold", justify="CT", no_clip=True)
    fig.savefig("fig_classfrac.png", dpi=300)
    fig.savefig("fig_classfrac.pdf")

def main():
    res = run_test(*load())
    for c in CLASSES:
        d = res[c]
        print("%-14s nb=%d obs=%.1f null=%.1f p=%.4f frac=%.3f"
              % (c, d["nb"], d["obs"], d["nulls"].mean(), d["p"], d["fobs"]))
    fig_coincidence(res)
    fig_classfrac(res)
    print("saved fig_coincidence.{png,pdf} and fig_classfrac.{png,pdf}")

if __name__ == "__main__":
    main()
