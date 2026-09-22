#!/usr/bin/env bash
set -e
AXIS="axes_full/1_Aleutian.txt"; GRID="@earth_relief_02m"

python3 - "$AXIS" << 'PY'
import math, json, sys
R=6371.0; DS=15.0; L=50.0; DU=5.0
AX=[tuple(map(float,ln.split()[:2])) for ln in open(sys.argv[1]) if ln.strip() and not ln.startswith("#")]
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
json.dump({"N":N,"DS":DS},open("meta.json","w"))
PY

gmt grdtrack samples.txt -G$GRID > sampled.txt

python3 << 'PY'
import json, math, statistics
m=json.load(open("meta.json")); N=m["N"]; DS=m["DS"]
best={}
for line in open("sampled.txt"):
    p=line.split()
    if len(p)<5 or p[4]=="NaN": continue
    idx=int(float(p[2])); z=float(p[4])
    if idx not in best or z<best[idx]: best[idx]=z
s=[i*DS for i in range(N) if i in best]; d=[best[i]/1000.0 for i in range(N) if i in best]
n=len(d); smax=s[-1]
def mv(a,w):
    h=w//2; return [sum(a[max(0,i-h):min(len(a),i+h+1)])/len(a[max(0,i-h):min(len(a),i+h+1)]) for i in range(len(a))]
def partition(dd,beta,LMIN=4):
    nn=len(dd); c1=[0.0]*(nn+1); c2=[0.0]*(nn+1)
    for k in range(nn): c1[k+1]=c1[k]+dd[k]; c2[k+1]=c2[k]+dd[k]**2
    def sc(a,b):
        m2=b-a; s1=c1[b]-c1[a]; s2=c2[b]-c2[a]; return s2-s1*s1/m2
    F=[None]*(nn+1); F[0]=-beta; prev=[0]*(nn+1)
    for t in range(1,nn+1):
        hi=t-LMIN; cs=range(0,hi+1) if hi>=0 else [0]
        bv=None; bs=0
        for a in cs:
            cc=F[a]+sc(a,t)+beta
            if bv is None or cc<bv: bv=cc; bs=a
        F[t]=bv; prev[t]=bs
    cps=[]; t=nn
    while t>0:
        a=prev[t]
        if a>0: cps.append(a)
        t=a
    return sorted(cps)
SMOOTHS=[3,5,7,9]; FACS=[1,1.5,2,2.7,3.5,4.5,6,7.5,9,11]; THR=0.6
COLORS={3:"gray55",5:"dodgerblue",7:"navy",9:"firebrick"}
allbounds=[]; lines_by_w={}
for w in SMOOTHS:
    dsm=mv(d,w); mu=sum(dsm)/n; sd=(sum((x-mu)**2 for x in dsm)/n)**0.5 or 1.0
    dd=[(x-mu)/sd for x in dsm]; sig2=0.5*statistics.pvariance([dd[k]-dd[k-1] for k in range(1,n)])
    pts=[]
    for f in FACS:
        beta=f*max(sig2,1e-6)*math.log(n); cps=partition(dd,beta); pts.append((f,len(cps)))
        for c in cps: allbounds.append(s[c])
    lines_by_w[w]=pts
    open(f"lineK_{w}.txt","w").write("\n".join(f"{f} {K}" for f,K in pts)+"\n")
ncomb=len(SMOOTHS)*len(FACS); binw=30.0; nb=int(smax//binw)+1
cnt=[0]*nb
for sb in allbounds: cnt[min(int(sb//binw),nb-1)]+=1
rec=[c/ncomb for c in cnt]
open("recur.txt","w").write("\n".join(f"{(k+0.5)*binw} {rec[k]}" for k in range(nb))+"\n")
stable=[]; k=0
while k<nb:
    if rec[k]>=THR:
        j=k
        while j+1<nb and rec[j+1]>=THR: j+=1
        bb=max(range(k,j+1),key=lambda b:rec[b]); stable.append((bb+0.5)*binw); k=j+1
    else: k+=1
open("stable.txt","w").write("\n".join(f"{x} 1.03" for x in stable)+"\n")
open("profileB.txt","w").write("\n".join(f"{s[k]} {mv(d,7)[k]}" for k in range(n))+"\n")
dmin=min(d)-0.3; dmax=max(d)+0.3
with open("legA.txt","w") as f:
    f.write("L 9p,Helvetica C Smoothing width (nodes)\nN 4\n")
    for w in SMOOTHS: f.write(f"S 0.3c - 1.0c - 1p,{COLORS[w]} 0.4c {w}\n")
kmax=max(K for w in SMOOTHS for _,K in lines_by_w[w])+2
open("info.txt","w").write(f"{smax:.0f} {min(FACS)} {max(FACS)} {kmax} {dmin:.2f} {dmax:.2f} {len(stable)} {ncomb} {THR}\n")
PY

read SMAX FMIN FMAX KMAX DMIN DMAX NST NC THR < info.txt
gmt begin fig_sensitivity pdf,png
  gmt set FONT_TITLE 11p,Helvetica FONT_ANNOT_PRIMARY 9p FONT_LABEL 10p MAP_FRAME_TYPE plain MAP_GRID_PEN_PRIMARY 0.2p,gray88
  gmt basemap -R0/$SMAX/0/1.1 -JX17c/5c -BWSne+t"(b) Boundary recurrence across the parameter grid (Aleutian)" \
      -Bxa500f100g500+l"Along-strike distance (km)" -Bya0.2f0.1g0.2+l"Recurrence (fraction of grid)"
  gmt plot recur.txt -Sb25u+b0 -Gsteelblue -W0.2p,gray30
  printf "0 %s\n%s %s\n" "$THR" "$SMAX" "$THR" | gmt plot -W0.4p,red,--
  gmt plot stable.txt -Si0.30c -Gred -W0.3p,black -N
  printf "1300 %s stable (>=%s)\n" "$THR" "$THR" | gmt text -F+f8p,Helvetica-Bold,red+jCB -D0/0.12c -Gwhite@25 -N
  gmt basemap -R0/$SMAX/$DMIN/$DMAX -JX17c/5c -BE -Bya2f1+l"Axial depth (km)"
  gmt plot profileB.txt -W0.8p,gray45
  gmt basemap -R$FMIN/$FMAX/0/$KMAX -JX17c/4.2c -BWSne+t"(a) Number of boundaries vs penalty, by smoothing width" \
      -Bxa2f1g1+l"Penalty factor (@~b@~ / BIC)" -Bya5f1g5+l"Boundaries K" -Y7.3c
  for w in 3 5 7 9; do
    case $w in 3) c=gray55;; 5) c=dodgerblue;; 7) c=navy;; 9) c=firebrick;; esac
    gmt plot lineK_${w}.txt -W0.5p,$c
    gmt plot lineK_${w}.txt -Sc0.08c -G$c
  done
  gmt legend legA.txt -DjBL+w9c+o0.3c/0.3c -F+gwhite@10+p0.4p,gray50
gmt end
echo "wrote fig_sensitivity.pdf and fig_sensitivity.png"
