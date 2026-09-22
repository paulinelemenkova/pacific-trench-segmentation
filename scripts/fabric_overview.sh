#!/usr/bin/env bash
set -e

REGION="115/295/-62/66"
PROJ="J200/18c"

python3 - << 'PY'
import glob, os, re
GROUP={"NW Pacific":[1,2,3,4,5],"W Pacific":[6,7,8,9,10],
       "SW Pacific":[11,12,13,14,15,16,17,18],"E Pacific":[19,20]}
GCOLOR={"NW Pacific":"firebrick","W Pacific":"darkorange","SW Pacific":"purple","E Pacific":"navy"}
tcolor={i:GCOLOR[g] for g,ids in GROUP.items() for i in ids}
tr={}
for p in glob.glob(os.path.join("axes_full","*.txt")):
    m=re.match(r"(\d+)_(.+)\.txt$", os.path.basename(p))
    if not m: continue
    tr[int(m.group(1))]=[ln.split()[:2] for ln in open(p) if ln.strip() and not ln.startswith("#")]
tr=dict(sorted(tr.items()))
with open("axes.gmt","w") as f, open("axes_plain.gmt","w") as g:
    for i,pts in tr.items():
        f.write(f"> -W1.6p,{tcolor[i]}\n"); g.write("> \n")
        for lo,la in pts:
            f.write(f"{lo} {la}\n"); g.write(f"{lo} {la}\n")

plateaus=[(160,-4,"Ontong Java Plateau"),(159,33,"Shatsky Rise"),(178,35,"Hess Rise"),
          (196,-10,"Manihiki Plateau"),(183,-41,"Hikurangi Plateau")]
ridge_above=[(171,44,"Emperor Smts"),(203,23,"Hawaiian Ridge"),(192,-33,"Louisville Ridge"),(154,19,"Marcus-Wake Smts")]
ridge_left=[(272,6,"Cocos Ridge"),(274,-2,"Carnegie Ridge"),(281,-18,"Nazca Ridge")]
fzs=[(228,40.5,"Mendocino FZ"),(231,18.5,"Clarion FZ"),(239,10.5,"Clipperton FZ")]
def dump(fn,rows,pts=False):
    with open(fn,"w") as f:
        for lo,la,t in rows: f.write(f"{lo} {la}\n" if pts else f"{lo} {la} {t}\n")
dump("plat.txt",plateaus); dump("plat_pts.txt",plateaus,True)
dump("ridge_above.txt",ridge_above); dump("ridge_left.txt",ridge_left)
dump("ridge_pts.txt",ridge_above+ridge_left,True); dump("fz.txt",fzs)
with open("legend.txt","w") as f:
    f.write("N 4\n")
    f.write("S 0.22c c 0.26c forestgreen 0.3p,black 0.55c Oceanic plateaus\n")
    f.write("S 0.22c a 0.32c magenta4 0.3p,black 0.55c Aseismic ridges / seamounts\n")
    f.write("G 0.05c\nN 2\n")
    f.write("S 0.25c - 0.7c - 2.2p,purple 0.8c Trench axes (by sector, see study-area map)\n")
    f.write("L 8p,Helvetica-Oblique L Fracture zones labelled in black italic\n")
PY

gmt begin fig_fabric_overview pdf,png
  gmt set FONT_TITLE 12p,Helvetica FONT_ANNOT_PRIMARY 10p \
          MAP_GRID_PEN_PRIMARY 0.5p,white MAP_FRAME_TYPE fancy
  gmt makecpt -Cvik -T-60/60/5 -H > faa.cpt
  gmt grdimage @earth_faa -R${REGION} -J${PROJ} -Cfaa.cpt -I+d \
      -Bxa30f10g30 -Bya20f10g20 -BWSne+t"Inherited seafloor fabric of the Pacific"
  gmt coast -Ggray82 -W0.3p,gray40 -Dl -A200
  gmt plot axes_plain.gmt -W2.6p,white
  gmt plot axes.gmt
  gmt plot plat_pts.txt  -Sc0.22c -Gforestgreen -W0.4p,white
  gmt plot ridge_pts.txt -Sa0.32c -Gmagenta4   -W0.4p,white
  gmt text plat.txt        -F+f8p,Helvetica-Bold,forestgreen+jCB -Gwhite@25 -W0.2p,gray -D0/0.30c -N
  gmt text ridge_above.txt -F+f8p,Helvetica-Bold,magenta4+jCB    -Gwhite@25 -W0.2p,gray -D0/0.30c -N
  gmt text ridge_left.txt  -F+f8p,Helvetica-Bold,magenta4+jRM    -Gwhite@25 -W0.2p,gray -D-0.30c/0 -N
  gmt text fz.txt          -F+f8p,Helvetica-Oblique,black+jCM    -Gwhite@25 -W0.2p,gray -N
  gmt colorbar -Cfaa.cpt -DjBR+w5c/0.32c+o0.7c/0.7c+v \
      -Bxa30f10+l"Free-air anomaly (mGal)" -F+gwhite@25+p0.4p,gray40
  gmt legend legend.txt -DJBC+jTC+w18c+o0/0.8c -F+gwhite@10+p0.5p,gray50
gmt end

echo "wrote fig_fabric_overview.pdf and fig_fabric_overview.png"
