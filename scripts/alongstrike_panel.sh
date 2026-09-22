#!/usr/bin/env bash
set -e
GRID="@earth_relief_04m"
DS=15; L=40; DU=5

python3 - "$GRID" << 'PY'
import math, json, glob, os, re, sys
R=6371.0; DS=15.0; L=40.0; DU=5.0
def hav(a,b):
    lo1,la1=map(math.radians,a); lo2,la2=map(math.radians,b)
    d=math.sin((la2-la1)/2)**2+math.cos(la1)*math.cos(la2)*math.sin((lo2-lo1)/2)**2
    return 2*R*math.asin(math.sqrt(d))
def bearing(a,b):
    lo1,la1=map(math.radians,a); lo2,la2=map(math.radians,b)
    y=math.sin(lo2-lo1)*math.cos(la2); x=math.cos(la1)*math.sin(la2)-math.sin(la1)*math.cos(la2)*math.cos(lo2-lo1)
    return math.atan2(y,x)
def dest(lon,lat,brg,dist):
    la1,lo1,dR=math.radians(lat),math.radians(lon),dist/R
    la2=math.asin(math.sin(la1)*math.cos(dR)+math.cos(la1)*math.sin(dR)*math.cos(brg))
    lo2=lo1+math.atan2(math.sin(brg)*math.sin(dR)*math.cos(la1),math.cos(dR)-math.sin(la1)*math.sin(la2))
    return (math.degrees(lo2)%360, math.degrees(la2))
TR={}
for p in glob.glob("axes_full/*.txt"):
    m=re.match(r"(\d+)_(.+)\.txt$", os.path.basename(p))
    if not m: continue
    TR[int(m.group(1))]=(m.group(2).replace("_"," "),
        [tuple(map(float,ln.split()[:2])) for ln in open(p) if ln.strip() and not ln.startswith("#")])
names={}
with open("samples.txt","w") as f:
    for tid in sorted(TR):
        name,AX=TR[tid]; names[tid]=name
        cum=[0.0]
        for i in range(1,len(AX)): cum.append(cum[-1]+hav(AX[i-1],AX[i]))
        dense=[]
        for s in [k*DS for k in range(int(cum[-1]//DS)+1)]:
            j=0
            while j<len(cum)-2 and cum[j+1]<s: j+=1
            t=(s-cum[j])/max(cum[j+1]-cum[j],1e-9)
            dense.append((AX[j][0]+t*(AX[j+1][0]-AX[j][0]),AX[j][1]+t*(AX[j+1][1]-AX[j][1])))
        n=len(dense)
        brg=[bearing(dense[max(i-1,0)],dense[min(i+1,n-1)]) for i in range(n)]
        offs=[-L+k*DU for k in range(int(2*L//DU)+1)]
        for i,(lon,lat) in enumerate(dense):
            for u in offs:
                b=brg[i]+math.pi/2 if u>=0 else brg[i]-math.pi/2
                plon,plat=dest(lon,lat,b,abs(u)); f.write(f"{plon} {plat} {tid} {i} {u}\n")
json.dump({"DS":DS,"names":names},open("meta.json","w"))
PY

gmt grdtrack samples.txt -G$GRID > sampled.txt

python3 << 'PY'
import json, math, statistics
m=json.load(open("meta.json")); DS=m["DS"]; NAMES={int(k):v for k,v in m["names"].items()}
best={}
for line in open("sampled.txt"):
    p=line.split()
    if len(p)<6 or p[5]=="NaN": continue
    tid=int(p[2]); nid=int(p[3]); z=float(p[5]); key=(tid,nid)
    if key not in best or z<best[key]: best[key]=z
def movavg(a,w=7):
    h=w//2; return [sum(a[max(0,i-h):min(len(a),i+h+1)])/len(a[max(0,i-h):min(len(a),i+h+1)]) for i in range(len(a))]
def partition(dd,beta,LMIN=5):
    n=len(dd); c1=[0.0]*(n+1); c2=[0.0]*(n+1)
    for k in range(n): c1[k+1]=c1[k]+dd[k]; c2[k+1]=c2[k]+dd[k]**2
    def sc(a,b):
        nn=b-a; s1=c1[b]-c1[a]; s2=c2[b]-c2[a]; return s2-s1*s1/nn
    F=[None]*(n+1); F[0]=-beta; prev=[0]*(n+1)
    for t in range(1,n+1):
        hi=t-LMIN; cands=range(0,hi+1) if hi>=0 else [0]
        bv=None; bs=0
        for a in cands:
            cc=F[a]+sc(a,t)+beta
            if bv is None or cc<bv: bv=cc; bs=a
        F[t]=bv; prev[t]=bs
    cps=[]; t=n
    while t>0:
        a=prev[t]
        if a>0: cps.append(a)
        t=a
    return sorted(cps)
meta=[]
for tid in sorted(NAMES):
    nids=sorted(k[1] for k in best if k[0]==tid)
    if not nids: continue
    s=[nid*DS for nid in nids]; d=[best[(tid,nid)]/1000.0 for nid in nids]
    d=movavg(d,7); n=len(d)
    mu=sum(d)/n; sd=(sum((x-mu)**2 for x in d)/n)**0.5 or 1.0; dd=[(x-mu)/sd for x in d]
    sig2=0.5*statistics.pvariance([dd[k]-dd[k-1] for k in range(1,n)]) if n>2 else 1.0
    beta=4.5*max(sig2,1e-6)*math.log(max(n,2))
    cps=partition(dd,beta)
    open(f"series_{tid}.txt","w").write("\n".join(f"{s[k]} {d[k]}" for k in range(n))+"\n")
    dmin=min(d)-0.3; dmax=max(d)+0.3
    open(f"bounds_{tid}.txt","w").write("\n".join(f"> \n{s[c]} {dmin}\n{s[c]} {dmax}" for c in cps)+"\n")
    meta.append((tid,NAMES[tid],max(s[-1],1),round(dmin,2),round(dmax,2)))
with open("panelmeta.txt","w") as f:
    for tid,name,smax,dmin,dmax in meta: f.write(f"{tid}|{name}|{smax:.0f}|{dmin}|{dmax}\n")
open("legend.txt","w").write("N 2\nS 0.7c - 1.1c - 0.8p,navy 1.4c Axial-depth profile d(s)\nS 0.7c - 1.1c - 0.4p,red,-- 1.4c Detected segment boundary\n")
PY

gmt begin fig_alongstrike_panel pdf,png
  gmt set FONT_HEADING 12p,Helvetica FONT_TITLE 9p,Helvetica FONT_ANNOT_PRIMARY 6p FONT_LABEL 7p MAP_FRAME_TYPE plain MAP_FRAME_PEN 0.6p MAP_TICK_LENGTH_PRIMARY 0.08c MAP_GRID_PEN_PRIMARY 0.15p,gray90
  gmt subplot begin 5x4 -Fs4.2c/2.85c -M0.35c/0.85c -T"Along-strike axial-depth profiles with detected segment boundaries"
  idx=0
  while IFS='|' read tid name smax dmin dmax; do
    gmt subplot set $idx
    gmt basemap -R0/$smax/$dmin/$dmax -JX4.2c/2.85c -Bxafg -Byafg -BWSne+t"$tid. $name"
    gmt plot bounds_${tid}.txt -W0.4p,red,--
    gmt plot series_${tid}.txt -W0.8p,navy
    idx=$((idx+1))
  done < panelmeta.txt
  gmt subplot end
  gmt legend legend.txt -DJBC+jTC+w12c+o0/0.7c -F+gwhite+p0.5p,gray50
gmt end
echo "wrote fig_alongstrike_panel.pdf and fig_alongstrike_panel.png"
