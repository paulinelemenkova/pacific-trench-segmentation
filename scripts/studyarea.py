#!/usr/bin/env python3
import glob
import os
import re

import numpy as np
import pygmt

AXES_DIR = "axes_full"
REGION = [115, 295, -62, 66]
PROJ = "J200/18c"

GROUP = {
    "NW Pacific": [1, 2, 3, 4, 5],
    "W Pacific":  [6, 7, 8, 9, 10],
    "SW Pacific": [11, 12, 13, 14, 15, 16, 17, 18],
    "E Pacific":  [19, 20],
}
GCOLOR = {"NW Pacific": "firebrick", "W Pacific": "darkorange",
          "SW Pacific": "purple", "E Pacific": "navy"}
tcolor = {i: GCOLOR[g] for g, ids in GROUP.items() for i in ids}

LIGHT_CPT = """# light bright bathy/topo, hinge at 0
-8000 65/105/225 -6000 79/129/230
-6000 79/129/230 -4000 100/149/237
-4000 100/149/237 -2000 141/182/240
-2000 141/182/240 -800 176/224/230
-800 176/224/230 -200 205/240/248
-200 205/240/248 0 224/255/255
0 200/232/180 300 224/220/165
300 224/220/165 1200 210/198/158
1200 210/198/158 3000 222/212/190
3000 222/212/190 6000 250/248/245
B 65/105/225
F 250/248/245
N 200/200/200
"""

def load_axes(folder):
    out = {}
    for path in glob.glob(os.path.join(folder, "*.txt")):
        m = re.match(r"(\d+)_(.+)\.txt$", os.path.basename(path))
        if not m:
            continue
        num, name = int(m.group(1)), m.group(2).replace("_", " ")
        out[num] = (name, np.atleast_2d(np.loadtxt(path)))
    return dict(sorted(out.items()))

def write_overlays(trenches):

    with open("axes.gmt", "w") as f, open("axes_plain.gmt", "w") as g:
        for i, (name, a) in trenches.items():
            f.write(f"> -W2.8p,{tcolor[i]}\n")
            g.write("> \n")
            for lo, la in a:
                f.write(f"{lo} {la}\n")
                g.write(f"{lo} {la}\n")

    leg = ["L 9p,Helvetica-Bold C Trench axes coloured by sector", "N 4"]
    for gname in ["NW Pacific", "W Pacific", "SW Pacific", "E Pacific"]:
        leg.append(f"S 0.2c r 0.5c/0.13c {GCOLOR[gname]} - 0.5c {gname}")
    leg += ["D 0.08c 0.6p", "N 4"]
    for i, (name, a) in trenches.items():
        leg.append(f"S 0.2c r 0.5c/0.13c {tcolor[i]} - 0.55c {i}. {name}")
    open("legend.txt", "w").write("\n".join(leg) + "\n")
    open("light_bathy.cpt", "w").write(LIGHT_CPT)

def main():
    trenches = load_axes(AXES_DIR)
    write_overlays(trenches)

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="fancy", FONT_TITLE="12p,Helvetica",
                 FONT_ANNOT_PRIMARY="10p", MAP_GRID_PEN_PRIMARY="0.5p,white")

    fig.grdimage(grid="@earth_relief_10m", region=REGION, projection=PROJ,
                 cmap="light_bathy.cpt", shading="+d",
                 frame=["xa30f10g30", "ya20f10g20",
                        "WSne+tPacific subduction trenches"])
    fig.coast(shorelines="0.2p,gray30", resolution="l")

    fig.plot(data="axes_plain.gmt", pen="4p,white")
    fig.plot(data="axes.gmt")

    fig.text(x=205, y=8, text="PACIFIC OCEAN", no_clip=True,
             font="18p,Helvetica-Bold,azure=0.5p,steelblue4")

    xs, ys, ids = [], [], []
    for i, (name, a) in trenches.items():
        mx, my = a[len(a) // 2]
        xs.append(mx); ys.append(my); ids.append(str(i))
    fig.text(x=xs, y=ys, text=ids, no_clip=True,
             font="6p,Helvetica-Bold,white", fill="black",
             pen="0.25p,black", clearance="0.03c/0.03c+tO")

    fig.legend(spec="legend.txt",
               position="JBC+jTC+w18c+o0/0.8c",
               box="+gwhite@10+p0.5p,gray50")

    fig.savefig("fig_studyarea.pdf")
    fig.savefig("fig_studyarea.png", dpi=300)
    print("wrote fig_studyarea.pdf and fig_studyarea.png")

if __name__ == "__main__":
    main()
