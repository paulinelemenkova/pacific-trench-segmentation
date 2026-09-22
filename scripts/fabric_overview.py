#!/usr/bin/env python3
import glob
import os
import re

import numpy as np
import pygmt

AXES_DIR = "axes_full"
REGION = [115, 295, -62, 66]
PROJ = "J200/18c"

GROUP = {"NW Pacific": [1, 2, 3, 4, 5], "W Pacific": [6, 7, 8, 9, 10],
         "SW Pacific": [11, 12, 13, 14, 15, 16, 17, 18], "E Pacific": [19, 20]}
GCOLOR = {"NW Pacific": "firebrick", "W Pacific": "darkorange",
          "SW Pacific": "purple", "E Pacific": "navy"}
tcolor = {i: GCOLOR[g] for g, ids in GROUP.items() for i in ids}

PLATEAUS = [(160, -4, "Ontong Java Plateau"), (159, 33, "Shatsky Rise"),
            (178, 35, "Hess Rise"), (196, -10, "Manihiki Plateau"),
            (183, -41, "Hikurangi Plateau")]
RIDGE_ABOVE = [(171, 44, "Emperor Smts"), (203, 23, "Hawaiian Ridge"),
               (192, -33, "Louisville Ridge"), (154, 19, "Marcus-Wake Smts")]
RIDGE_LEFT = [(272, 6, "Cocos Ridge"), (274, -2, "Carnegie Ridge"),
              (281, -18, "Nazca Ridge")]
FZS = [(228, 40.5, "Mendocino FZ"), (231, 18.5, "Clarion FZ"),
       (239, 10.5, "Clipperton FZ")]

def load_axes(folder):
    out = {}
    for path in glob.glob(os.path.join(folder, "*.txt")):
        m = re.match(r"(\d+)_(.+)\.txt$", os.path.basename(path))
        if not m:
            continue
        out[int(m.group(1))] = np.atleast_2d(np.loadtxt(path))
    return dict(sorted(out.items()))

def write_axes(trenches):
    with open("axes.gmt", "w") as f, open("axes_plain.gmt", "w") as g:
        for i, a in trenches.items():
            f.write(f"> -W1.6p,{tcolor[i]}\n")
            g.write("> \n")
            for lo, la in a:
                f.write(f"{lo} {la}\n")
                g.write(f"{lo} {la}\n")

def labels(fig, rows, color, justify, offset):
    xs = [r[0] for r in rows]; ys = [r[1] for r in rows]; ts = [r[2] for r in rows]
    fig.text(x=xs, y=ys, text=ts, font=f"8p,Helvetica-Bold,{color}",
             justify=justify, offset=offset, fill="white@25",
             pen="0.2p,gray", no_clip=True)

def main():
    trenches = load_axes(AXES_DIR)
    write_axes(trenches)

    leg = ["N 4",
           "S 0.22c c 0.26c forestgreen 0.3p,black 0.55c Oceanic plateaus",
           "S 0.22c a 0.32c magenta4 0.3p,black 0.55c Aseismic ridges / seamounts",
           "G 0.05c", "N 2",
           "S 0.25c - 0.7c - 2.2p,purple 0.8c Trench axes (by sector, see study-area map)",
           "L 8p,Helvetica-Oblique L Fracture zones labelled in black italic"]
    open("legend.txt", "w").write("\n".join(leg) + "\n")

    pygmt.makecpt(cmap="vik", series=[-60, 60, 5], output="faa.cpt")

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="fancy", FONT_TITLE="12p,Helvetica",
                 FONT_ANNOT_PRIMARY="10p", MAP_GRID_PEN_PRIMARY="0.5p,white")

    fig.grdimage(grid="@earth_faa", region=REGION, projection=PROJ,
                 cmap="faa.cpt", shading="+d",
                 frame=["xa30f10g30", "ya20f10g20",
                        "WSne+tInherited seafloor fabric of the Pacific"])
    fig.coast(land="gray82", shorelines="0.3p,gray40", resolution="l", area_thresh=200)

    fig.plot(data="axes_plain.gmt", pen="2.6p,white")
    fig.plot(data="axes.gmt")

    fig.plot(x=[r[0] for r in PLATEAUS], y=[r[1] for r in PLATEAUS],
             style="c0.22c", fill="forestgreen", pen="0.4p,white")
    rp = RIDGE_ABOVE + RIDGE_LEFT
    fig.plot(x=[r[0] for r in rp], y=[r[1] for r in rp],
             style="a0.32c", fill="magenta4", pen="0.4p,white")

    labels(fig, PLATEAUS, "forestgreen", "CB", "0/0.30c")
    labels(fig, RIDGE_ABOVE, "magenta4", "CB", "0/0.30c")
    labels(fig, RIDGE_LEFT, "magenta4", "RM", "-0.30c/0")
    xs = [r[0] for r in FZS]; ys = [r[1] for r in FZS]; ts = [r[2] for r in FZS]
    fig.text(x=xs, y=ys, text=ts, font="8p,Helvetica-Oblique,black",
             justify="CM", fill="white@25", pen="0.2p,gray", no_clip=True)

    fig.colorbar(cmap="faa.cpt",
                 position="jBR+w5c/0.32c+o0.7c/0.7c+v",
                 frame=["xa30f10", "+lFree-air anomaly (mGal)"],
                 box="+gwhite@25+p0.4p,gray40")
    fig.legend(spec="legend.txt", position="JBC+jTC+w18c+o0/0.8c",
               box="+gwhite@10+p0.5p,gray50")

    fig.savefig("fig_fabric_overview.pdf")
    fig.savefig("fig_fabric_overview.png", dpi=300)
    print("wrote fig_fabric_overview.pdf and fig_fabric_overview.png")

if __name__ == "__main__":
    main()
