#!/usr/bin/env bash
set -e
AXIS="axes_full/1_Aleutian.txt"
GRID="@earth_relief_02m"

python3 - "$AXIS" << 'PY'
import math, json, sys
R=6371.0; DS=12.0; L=50.0; DU=4.0
AX=[tuple(map(float,ln.split()[:2])) for ln in open(sys.argv[1]) if ln.strip() and not ln.startswith("#")]
def hav(a,b):
    lo1,la1=map(math.radians,a); lo2,la2=map(math.radians,b)
    d=math.sin((la2-la1)/2)**2+math.cos(la1)*math.cos(la2)*math.sin((lo2-lo1)/2)**2
    return 2*R*math.asin(math.sqrt(d))
def bearing(a,b):
    lo1,la1=map(math.radians,a); lo2,la2=map(math.radians,b)
    y=math.sin(lo2-lo1)*math.cos(la2)
    x=math.cos(la1)*math.sin(la2)-math.sin(la1)*math.cos(la2)*math.cos(lo2-lo1)
    return math.atan2(y,x)
def dest(lon,lat,brg,dist):
    la1,lo1,dR=math.radians(lat),math.radians(lon),dist/R
    la2=math.asin(math.sin(la1)*math.cos(dR)+math.cos(la1)*math.sin(dR)*math.cos(brg))
    lo2=lo1+math.atan2(math.sin(brg)*math.sin(dR)*math.cos(la1),math.cos(dR)-math.sin(la1)*math.sin(la2))
    return (math.degrees(lo2)%360, math.degrees(la2))
cum=[0.0]
for i in range(1,len(AX)): cum.append(cum[-1]+hav(AX[i-1],AX[i]))
dense=[]
for s in [k*DS for k in range(int(cum[-1]//DS)+1)]:
    j=0
    while j<len(cum)-2 and cum[j+1]<s: j+=1
    t=(s-cum[j])/max(cum[j+1]-cum[j],1e-9)
    dense.append((AX[j][0]+t*(AX[j+1][0]-AX[j][0]),AX[j][1]+t*(AX[j+1][1]-AX[j][1])))
N=len(dense)
brg=[bearing(dense[max(i-1,0)],dense[min(i+1,N-1)]) for i in range(N)]
offs=[-L+k*DU for k in range(int(2*L//DU)+1)]
with open("samples.txt","w") as f:
    for i,(lon,lat) in enumerate(dense):
        for u in offs:
            b=brg[i]+math.pi/2 if u>=0 else brg[i]-math.pi/2
            plon,plat=dest(lon,lat,b,abs(u)); f.write(f"{plon} {plat} {i} {u}\n")
open("axis_dense.txt","w").write("\n".join(f"{lo} {la}" for lo,la in dense)+"\n")
json.dump({"N":N,"DS":DS,"noff":len(offs)},open("meta.json","w"))
PY

gmt grdtrack samples.txt -G$GRID > sampled.txt

python3 << 'PY'
import json
m=json.load(open("meta.json")); N=m["N"]; DS=m["DS"]; SMOOTH=9
best={}; prof={}
for line in open("sampled.txt"):
    p=line.split()
    if len(p)<5 or p[4]=="NaN": continue
    idx=int(float(p[2])); z=float(p[4])
    prof.setdefault(idx,[]).append((float(p[0]),float(p[1]),z))
    if idx not in best or z<best[idx][0]: best[idx]=(z,float(p[0]),float(p[1]))
s=[];d=[];pk=[]
for i in range(N):
    if i in best: s.append(i*DS); d.append(best[i][0]/1000.0); pk.append((best[i][1],best[i][2]))
def mv(a,w):
    h=w//2; return [sum(a[max(0,i-h):min(len(a),i+h+1)])/len(a[max(0,i-h):min(len(a),i+h+1)]) for i in range(len(a))]
ds=mv(d,SMOOTH)
open("series_raw.txt","w").write("\n".join(f"{x} {y}" for x,y in zip(s,d))+"\n")
open("series_smooth.txt","w").write("\n".join(f"{x} {y}" for x,y in zip(s,ds))+"\n")
open("picks.txt","w").write("\n".join(f"{a} {b}" for a,b in pk)+"\n")
step=max(1,N//26)
with open("fan.txt","w") as f:
    for i in range(0,N,step):
        pl=prof.get(i)
        if pl: (lo1,la1,_),(lo2,la2,_)=pl[0],pl[-1]; f.write(f"> \n{lo1} {la1}\n{lo2} {la2}\n")
open("legend.txt","w").write("N 3\nS 0.4c - 0.9c - 1.4p,gold1,-- 1.1c Digitised axis\nS 0.4c - 0.9c - 1.8p,red 1.1c Deepest-point locus (Eq. 1)\nS 0.4c - 0.9c - 0.6p,chartreuse1 1.1c Cross-profiles (fan)\n")
open("axinfo.txt","w").write(f"{s[-1]:.0f} {round(min(d)-0.3,2)} {round(max(d)+0.3,2)}\n")
PY

read SMAX ZLO ZHI < axinfo.txt
gmt begin fig_axis_extraction pdf,png
  gmt set FONT_TITLE 11p,Helvetica FONT_ANNOT_PRIMARY 9p FONT_LABEL 10p MAP_FRAME_TYPE plain MAP_GRID_PEN_PRIMARY 0.25p,gray80
  gmt basemap -R0/$SMAX/$ZLO/$ZHI -JX18c/5c -BWSne+t"(b) Along-strike axial-depth series d(s)" \
      -Bxa500f100g500+l"Along-strike distance (km)" -Bya1f0.5g1+l"Axial depth (km)"
  gmt plot series_raw.txt -W0.6p,gray70
  gmt plot series_smooth.txt -W1p,navy
  gmt makecpt -Cgeo -T-8000/5000/100
  gmt set MAP_GRID_PEN_PRIMARY 0.2p,white
  gmt grdimage $GRID -R162/217/49.5/61 -JM18c -I+d -Y7.0c \
      -Bxa10f5g10 -Bya4f2g4 -BWsNe+t"(a) Cross-profile sampling and deepest-point picks"
  gmt coast -Ggray82 -W0.3p,gray45 -Di
  gmt plot fan.txt -W0.5p,chartreuse1
  gmt plot axis_dense.txt -W1.4p,gold1,--
  gmt plot picks.txt -W1.8p,red
  gmt legend legend.txt -DJBC+jTC+w15c+o0/0.35c -F+gwhite@10+p0.5p,gray50
gmt end
echo "wrote fig_axis_extraction.pdf and fig_axis_extraction.png"
