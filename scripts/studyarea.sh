#!/usr/bin/env bash
set -e

REGION="115/295/-62/66"
PROJ="J200/18c"

cat > light_bathy.cpt << 'CPT'
# light bright bathy/topo, hinge at 0
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
CPT

python3 - << 'PY'
import glob, os, re
GROUP = {"NW Pacific":[1,2,3,4,5], "W Pacific":[6,7,8,9,10],
         "SW Pacific":[11,12,13,14,15,16,17,18], "E Pacific":[19,20]}
GCOLOR = {"NW Pacific":"firebrick","W Pacific":"darkorange",
          "SW Pacific":"purple","E Pacific":"navy"}
tcolor = {i:GCOLOR[g] for g,ids in GROUP.items() for i in ids}
tr = {}
for p in glob.glob(os.path.join("axes_full","*.txt")):
    m = re.match(r"(\d+)_(.+)\.txt$", os.path.basename(p))
    if not m: continue
    num, name = int(m.group(1)), m.group(2).replace("_"," ")
    pts = [ln.split()[:2] for ln in open(p) if ln.strip() and not ln.startswith("#")]
    tr[num] = (name, pts)
tr = dict(sorted(tr.items()))
with open("axes.gmt","w") as f, open("axes_plain.gmt","w") as g:
    for i,(name,pts) in tr.items():
        f.write(f"> -W2.8p,{tcolor[i]}\n"); g.write("> \n")
        for lo,la in pts:
            f.write(f"{lo} {la}\n"); g.write(f"{lo} {la}\n")
with open("labels.txt","w") as f:
    for i,(name,pts) in tr.items():
        mx,my = pts[len(pts)//2]; f.write(f"{mx} {my} {i}\n")
with open("legend.txt","w") as f:
    f.write("L 9p,Helvetica-Bold C Trench axes coloured by sector\nN 4\n")
    for gname in ["NW Pacific","W Pacific","SW Pacific","E Pacific"]:
        f.write(f"S 0.2c r 0.5c/0.13c {GCOLOR[gname]} - 0.5c {gname}\n")
    f.write("D 0.08c 0.6p\nN 4\n")
    for i,(name,pts) in tr.items():
        f.write(f"S 0.2c r 0.5c/0.13c {tcolor[i]} - 0.55c {i}. {name}\n")
PY

gmt begin fig_studyarea pdf,png
  gmt set FONT_TITLE 12p,Helvetica FONT_ANNOT_PRIMARY 10p \
          MAP_GRID_PEN_PRIMARY 0.5p,white MAP_FRAME_TYPE fancy
  gmt grdimage @earth_relief_10m -R${REGION} -J${PROJ} -Clight_bathy.cpt -I+d \
      -Bxa30f10g30 -Bya20f10g20 -BWSne+t"Pacific subduction trenches"
  gmt coast -W0.2p,gray30 -Dl
  gmt plot axes_plain.gmt -W4.0p,white
  gmt plot axes.gmt
  echo "205 8 PACIFIC OCEAN" | \
      gmt text -F+f18p,Helvetica-Bold,azure=0.5p,steelblue4+a0 -N
  gmt text labels.txt -F+f6p,Helvetica-Bold,white -Gblack -W0.25p,black \
      -C0.03c/0.03c -N
  gmt legend legend.txt -DJBC+jTC+w18c+o0/0.8c -F+gwhite@10+p0.5p,gray50
gmt end

echo "wrote fig_studyarea.pdf and fig_studyarea.png"
