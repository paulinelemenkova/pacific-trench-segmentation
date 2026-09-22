#!/usr/bin/env bash
set -e

python3 << 'PY'
import csv, glob, os, re
DTOL=15.0; COL={"fracture_zone":"fz","plateau":"pl","ridge":"rg"}
# axes -> multisegment file
with open("axes.gmt","w") as f:
    _paths=[p for p in glob.glob("axes_full/*.txt") if re.match(r"\d+_",os.path.basename(p))]
    for p in sorted(_paths,
                    key=lambda s:int(re.match(r"(\d+)_",os.path.basename(s)).group(1))):
        f.write("> \n")
        for ln in open(p):
            if ln.strip() and not ln.startswith("#"):
                x,y=ln.split()[:2]; f.write(f"{float(x)%360} {y}\n")
bnd={}
for r in csv.DictReader(open("boundary_catalogue.csv")):
    bnd.setdefault(int(r["trench_id"]),[]).append((float(r["s_km"]),float(r["lon"]),float(r["lat"])))
con={}
for r in csv.DictReader(open("fabric_intersections.csv")):
    con.setdefault(int(r["trench_id"]),{}).setdefault(r["class"],[]).append(float(r["s_km"]))
pts={"bnone":[],"bfz":[],"bpl":[],"brg":[]}
for tid,bl in bnd.items():
    cc=con.get(tid,{})
    for s,lo,la in bl:
        best=bc=None
        for c in COL:
            if cc.get(c):
                d=min(abs(s-f) for f in cc[c])
                if best is None or d<best: best,bc=d,c
        key=("b"+COL[bc]) if (bc and best<=DTOL) else "bnone"
        pts[key].append(f"{lo} {la}")
for k,v in pts.items(): open(f"{k}.xy","w").write("\n".join(v)+"\n")
fab={c:[] for c in COL}
for r in csv.DictReader(open("fabric_intersections.csv")):
    fab[r["class"]].append(f'{r["lon"]} {r["lat"]}')
for c,ab in COL.items(): open(f"f{ab}.xy","w").write("\n".join(fab[c])+"\n")
open("gray.cpt","w").write("-8000 105 -3000 150\n-3000 150 -1000 185\n-1000 185 0 205\n"
    "0 205 1500 225\n1500 225 6000 245\nB 105\nF 245\nN 200\n")
open("legend.txt","w").write(
 "N 2\n"
 "S 0.2c c 0.12c gold 0.2p,black 0.55c Detected segment boundary\n"
 "S 0.2c d 0.17c chartreuse 0.35p,black 0.55c Fabric contact\n"
 "D 0.06c 0.5p\n"
 "L 8p,Helvetica-Bold C Boundary coinciding with fabric (@~\\243@~15 km):\n"
 "N 3\n"
 "S 0.2c c 0.22c cyan 0.4p,black 0.7c Fracture zone\n"
 "S 0.2c c 0.22c red1 0.4p,black 0.7c Oceanic plateau\n"
 "S 0.2c c 0.22c magenta 0.4p,black 0.7c Seamount / ridge\n")
# named inherited-fabric features to annotate: lon lat justify dx/dy  label
open("fab_labels.txt","w").write("\n".join([
 "167   53    BR -0.1c/0.55c Emperor Smt. Chain / Meiji Smt.",
 "188   51.5  BL 0.1c/0.55c Amlia FZ",
 "170   52.5  TC 0c/-0.55c Stalemate FZ",
 "158   33    ML 0.55c/0c Shatsky Rise",
 "152   20    ML 0.55c/0c Marcus-Wake Smts.",
 "138   8     TL 0.4c/-0.4c Caroline Ridge",
 "156   -5    ML 0.55c/0.15c Ontong Java Plateau",
 "167   -15   ML 0.55c/0c D'Entrecasteaux Ridge",
 "186   -26   ML 0.5c/-0.1c Louisville Smt. Chain",
 "179   -40   TC 0c/-0.55c Hikurangi Plateau",
 "276   8.5   MR -0.5c/0.15c Cocos Ridge",
 "279   -1    MR -0.5c/-0.1c Carnegie Ridge",
 "284   -15.5 MR -0.5c/0c Nazca Ridge",
 "281   -33   MR -0.5c/-0.1c Juan Fernandez Ridge",
])+"\n")
PY

gmt begin fig_boundary_map png,pdf
  gmt set FONT_TITLE 13p,Helvetica-Bold FONT_ANNOT_PRIMARY 9p MAP_FRAME_TYPE fancy MAP_GRID_PEN_PRIMARY 0.4p,white
  gmt grdimage @earth_relief_10m -R115/295/-62/66 -JJ200/18c -Cgray.cpt -I+d \
      -Bxa30f10g30 -Bya20f10g20 -BWSne+t"Detected segment boundaries and coinciding subducting fabric"
  gmt coast -Gwhite -W0.2p,gray55 -Dl -A200
  gmt plot axes.gmt -W1.0p,gray35
  [ -s ffz.xy ] && gmt plot ffz.xy -Sd0.17c -Gchartreuse -W0.35p,black
  [ -s fpl.xy ] && gmt plot fpl.xy -Sd0.17c -Gchartreuse -W0.35p,black
  [ -s frg.xy ] && gmt plot frg.xy -Sd0.17c -Gchartreuse -W0.35p,black
  [ -s bnone.xy ] && gmt plot bnone.xy -Sc0.06c -Ggold -W0.15p,black
  [ -s bfz.xy ] && gmt plot bfz.xy -Sc0.22c -Gcyan    -W0.5p,black
  [ -s bpl.xy ] && gmt plot bpl.xy -Sc0.22c -Gred1    -W0.5p,black
  [ -s brg.xy ] && gmt plot brg.xy -Sc0.22c -Gmagenta -W0.5p,black
  while read -r lon lat just off label; do
    [ -z "$lon" ] && continue
    echo "$lon $lat $label" | gmt text -F+f8p,Helvetica-Bold,black+j$just \
        -Gwhite -C0.03c -D$off+v0.4p,black
  done < fab_labels.txt
  gmt legend legend.txt -DJBC+jTC+w17c+o0/0.7c -F+gwhite@10+p0.5p,gray50
gmt end
echo "saved fig_boundary_map.png / .pdf"
