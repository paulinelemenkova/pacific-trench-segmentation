#!/usr/bin/env python3
import csv, glob, os, re
import pygmt

DTOL   = 15.0
REGION = [115, 295, -62, 66]
PROJ   = "J200/18c"
AXDIR  = "axes_full"
COL    = {"fracture_zone": "cyan", "plateau": "red1", "ridge": "magenta"}
AB     = {"fracture_zone": "fz", "plateau": "pl", "ridge": "rg"}
C_BOUNDARY = "gold"
C_CONTACT  = "chartreuse"

FABRIC_LABELS = [
    (167,  53.0, "Emperor Smt. Chain / Meiji Smt.", "BR", "-0.1c/0.55c"),
    (188,  51.5, "Amlia FZ",                          "BL", "0.1c/0.55c"),
    (170,  52.5, "Stalemate FZ",                      "TC", "0c/-0.55c"),
    (158,  33.0, "Shatsky Rise",                      "ML", "0.55c/0c"),
    (152,  20.0, "Marcus-Wake Smts.",                 "ML", "0.55c/0c"),
    (138,   8.0, "Caroline Ridge",                    "TL", "0.4c/-0.4c"),
    (156,  -5.0, "Ontong Java Plateau",               "ML", "0.55c/0.15c"),
    (167, -15.0, "D'Entrecasteaux Ridge",             "ML", "0.55c/0c"),
    (186, -26.0, "Louisville Smt. Chain",             "ML", "0.5c/-0.1c"),
    (179, -40.0, "Hikurangi Plateau",                 "TC", "0c/-0.55c"),
    (276,   8.5, "Cocos Ridge",                       "MR", "-0.5c/0.15c"),
    (279,  -1.0, "Carnegie Ridge",                    "MR", "-0.5c/-0.1c"),
    (284, -15.5, "Nazca Ridge",                       "MR", "-0.5c/0c"),
    (281, -33.0, "Juan Fernandez Ridge",              "MR", "-0.5c/-0.1c"),
]

def load_axes():
    segs = []

    paths = [p for p in glob.glob(f"{AXDIR}/*.txt")
             if re.match(r"\d+_", os.path.basename(p))]
    for p in sorted(paths,
                    key=lambda s: int(re.match(r"(\d+)_", os.path.basename(s)).group(1))):
        pts = [ln.split()[:2] for ln in open(p) if ln.strip() and not ln.startswith("#")]
        segs.append([(float(x) % 360, float(y)) for x, y in pts])
    return segs

def load():
    bnd = {}
    for r in csv.DictReader(open("boundary_catalogue.csv")):
        bnd.setdefault(int(r["trench_id"]), []).append(
            (float(r["s_km"]), float(r["lon"]), float(r["lat"])))
    con = {}
    for r in csv.DictReader(open("fabric_intersections.csv")):
        con.setdefault(int(r["trench_id"]), {}).setdefault(r["class"], []).append(float(r["s_km"]))
    return bnd, con

def classify(bnd, con):
    pts = {"bnone": [], "bfz": [], "bpl": [], "brg": []}
    for tid, bl in bnd.items():
        cc = con.get(tid, {})
        for s, lo, la in bl:
            best = bc = None
            for c in COL:
                if cc.get(c):
                    d = min(abs(s - f) for f in cc[c])
                    if best is None or d < best:
                        best, bc = d, c
            key = ("b" + AB[bc]) if (bc and best <= DTOL) else "bnone"
            pts[key].append((lo, la))
    return pts

def main():
    bnd, con = load()
    pts = classify(bnd, con)
    fab = {c: [] for c in COL}
    for r in csv.DictReader(open("fabric_intersections.csv")):
        fab[r["class"]].append((float(r["lon"]), float(r["lat"])))

    with open("_gray.cpt", "w") as f:
        f.write("-8000 105 -3000 150\n-3000 150 -1000 185\n-1000 185 0 205\n"
                "0 205 1500 225\n1500 225 6000 245\nB 105\nF 245\nN 200\n")
    with open("_legend.txt", "w") as f:
        f.write("N 2\n"
                "S 0.2c c 0.12c gold 0.2p,black 0.55c Detected segment boundary\n"
                "S 0.2c d 0.17c chartreuse 0.35p,black 0.55c Fabric contact\n"
                "D 0.06c 0.5p\n"
                "L 8p,Helvetica-Bold C Boundary coinciding with fabric (@~\\243@~15 km):\n"
                "N 3\n"
                "S 0.2c c 0.22c cyan 0.4p,black 0.7c Fracture zone\n"
                "S 0.2c c 0.22c red1 0.4p,black 0.7c Oceanic plateau\n"
                "S 0.2c c 0.22c magenta 0.4p,black 0.7c Seamount / ridge\n")

    fig = pygmt.Figure()
    pygmt.config(FONT_TITLE="13p,Helvetica-Bold", FONT_ANNOT_PRIMARY="9p",
                 MAP_FRAME_TYPE="fancy", MAP_GRID_PEN_PRIMARY="0.4p,white")
    fig.grdimage("@earth_relief_10m", region=REGION, projection=PROJ, cmap="_gray.cpt",
                 shading="+d", frame=["WSne+tDetected segment boundaries and coinciding "
                 "subducting fabric", "xa30f10g30", "ya20f10g20"])
    fig.coast(land="white", shorelines="0.2p,gray55", resolution="l", area_thresh=200)
    for s in load_axes():
        fig.plot(x=[p[0] for p in s], y=[p[1] for p in s], pen="1.0p,gray35")
    for c in COL:
        if fab[c]:
            fig.plot(x=[p[0] for p in fab[c]], y=[p[1] for p in fab[c]],
                     style="d0.17c", fill=C_CONTACT, pen="0.35p,black")
    if pts["bnone"]:
        fig.plot(x=[p[0] for p in pts["bnone"]], y=[p[1] for p in pts["bnone"]],
                 style="c0.06c", fill=C_BOUNDARY, pen="0.15p,black")
    for c, key in [("fracture_zone", "bfz"), ("plateau", "bpl"), ("ridge", "brg")]:
        if pts[key]:
            fig.plot(x=[p[0] for p in pts[key]], y=[p[1] for p in pts[key]],
                     style="c0.22c", fill=COL[c], pen="0.5p,black")

    for lon, lat, txt, just, off in FABRIC_LABELS:
        fig.text(x=lon, y=lat, text=txt, font="8p,Helvetica-Bold,black",
                 justify=just, fill="white", clearance="0.03c",
                 offset=f"{off}+v0.4p,black")
    fig.legend(spec="_legend.txt", position="JBC+jTC+w17c+o0/0.7c",
               box="+gwhite@10+p0.5p,gray50")
    fig.savefig("fig_boundary_map.png", dpi=300)
    fig.savefig("fig_boundary_map.pdf")
    print("saved fig_boundary_map.png / .pdf  |",
          {k: len(v) for k, v in pts.items()}, "contacts", {c: len(fab[c]) for c in COL})

if __name__ == "__main__":
    main()
