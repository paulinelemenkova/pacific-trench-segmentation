#!/usr/bin/env bash
set -e
IN=${IN:-.}
OUT=${OUT:-.}
REGION=115/288/-58/63
PROJ=JJ201.5/17c

cat > plat.txt <<EOF
158.5 33.0 CB Shatsky Rise
177.0 36.5 CB Hess Rise
159.0 -3.5 CB Ontong Java Plateau
197.0 -10.0 CB Manihiki Plateau
181.5 -41.0 LM Hikurangi Plateau
EOF
cat > ridge.txt <<EOF
170.5 44.5 RM Emperor Smts
198.0 22.5 LM Hawaiian Ridge
153.0 19.0 CB Marcus-Wake Smts
192.0 -33.5 LM Louisville Ridge
273.0 6.0 RM Cocos Ridge
274.5 -1.0 RM Carnegie Ridge
280.5 -18.0 RM Nazca Ridge
EOF
cat > fz.txt <<EOF
222.0 40.6 Mendocino FZ
233.0 18.8 Clarion FZ
238.0 10.3 Clipperton FZ
EOF

cat > legend.txt <<EOF
N 3
S 0.25c c 0.24c forestgreen 0.4p,white 0.55c Oceanic plateau
S 0.25c a 0.32c magenta4 0.4p,white 0.55c Aseismic ridge / seamount chain
S 0.25c - 0.001c white 0p,white 0.1c @%2%Fracture zones labelled in italic@%%
S 0.3c - 0.6c - 2.2p,yellow1 0.7c Trench axis, NW sector
S 0.3c - 0.6c - 2.2p,firebrick1 0.7c Trench axis, W sector
S 0.3c - 0.6c - 2.2p,chartreuse1 0.7c Trench axis, SW sector
S 0.3c - 0.6c - 2.2p,magenta 0.7c Trench axis, E sector
S 0.3c - 0.6c - 1.8p,black,6_3:0 0.7c Excluded trough (not analysed)
EOF

gmt begin "$OUT/fig_gravity_fabric" png,pdf E300
  gmt set FONT_ANNOT_PRIMARY 9p,Helvetica FONT_LABEL 9p,Helvetica \
          MAP_FRAME_TYPE plain MAP_GRID_PEN_PRIMARY 0.4p,white@30 \
          FORMAT_GEO_MAP ddd:mmF MAP_FRAME_PEN 1p,black
  H=$(gmt mapproject -R$REGION -$PROJ -Wh)
  echo MAPH $H
  gmt grdcut @earth_faa_05m -R$REGION -Gfaa.nc
  gmt makecpt -Cvik -T-60/60 -D -H > faa.cpt
  gmt grdimage faa.nc -R$REGION -$PROJ -Cfaa.cpt -I+a315+ne0.6 \
      -Bxa30f10g30 -Bya20f10g20 -BWSne
  gmt coast -Ggray82 -W0.25p,gray40 -Dl -A500
  for s in NW W SW E; do gmt plot "$IN/tr_$s.txt" -W3.6p,black; done
  for s in SWx Ex; do gmt plot "$IN/tr_$s.txt" -W3.4p,white; done
  gmt plot "$IN/tr_NW.txt" -W2.2p,yellow1
  gmt plot "$IN/tr_W.txt"  -W2.2p,firebrick1
  gmt plot "$IN/tr_SW.txt" -W2.2p,chartreuse1
  gmt plot "$IN/tr_E.txt"  -W2.2p,magenta
  gmt plot "$IN/tr_SWx.txt" -W1.8p,black,6_3:0
  gmt plot "$IN/tr_Ex.txt"  -W1.8p,black,6_3:0
  gmt plot plat.txt  -Sc0.24c -Gforestgreen -W0.4p,white
  gmt plot ridge.txt -Sa0.34c -Gmagenta4 -W0.4p,white
  gmt text plat.txt  -F+f9p,Helvetica-Bold,forestgreen+j -Gwhite@25 -C0.06c -Dj0.26c -N
  gmt text ridge.txt -F+f9p,Helvetica-Bold,magenta4+j -Gwhite@25 -C0.06c -Dj0.26c -N
  gmt text fz.txt    -F+f9p,Helvetica-Oblique,black+jCM -Gwhite@25 -C0.06c -N
  gmt colorbar -Cfaa.cpt -DJMR+w${H}c/0.3c+o0.35c/0+e -Bxa30f10 -By+l"mGal" \
      -Bx+l"Free-air anomaly"
  gmt legend legend.txt -DJBC+jTC+w17c+o0/0.75c -F+gwhite+p0.5p,gray50
gmt end
