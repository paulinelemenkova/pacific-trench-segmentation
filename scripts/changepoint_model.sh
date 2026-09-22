#!/usr/bin/env bash
set -e
AXIS="axes_full/1_Aleutian.txt"
GRID="@earth_relief_02m"

python3 - "$AXIS" << 'PY'
import math, json, sys
R=6371.0; DS=15.0; L=150.0; DU=3.0
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
prof={}
for line in open("sampled.txt"):
    p=line.split()
    if len(p)<5 or p[4]=="NaN": continue
    idx=int(float(p[2])); u=float(p[3]); z=float(p[4]); prof.setdefault(idx,[]).append((u,z))
s=[];d=[];W=[];A=[]
for i in range(N):
    pr=prof.get(i)
    if not pr or len(pr)<10: continue
    pr.sort()
    cand=[(z,k) for k,(u,z) in enumerate(pr) if abs(u)<=50]
    if not cand: continue
    zf,kf=min(cand); uf=pr[kf][0]; thr=zf+2000.0
    kL=kf
    while kL>0 and pr[kL][1]<thr: kL-=1
    kR=kf
    while kR<len(pr)-1 and pr[kR][1]<thr: kR+=1
    uL,zL=pr[kL]; uR,zR=pr[kR]
    gL=(zL-zf)/max(uf-uL,1e-6); gR=(zR-zf)/max(uR-uf,1e-6)
    s.append(i*DS); d.append(zf/1000.0); W.append(uR-uL); A.append((gL-gR)/max(gL+gR,1e-6))
n=len(s)
def zsc(a):
    mu=sum(a)/len(a); sd=(sum((x-mu)**2 for x in a)/len(a))**0.5 or 1.0
    return [(x-mu)/sd for x in a]
dd,WW,AA=zsc(d),zsc(W),zsc(A); Y=list(zip(dd,WW,AA)); p=3
cs1=[[0.0]*(n+1) for _ in range(p)]; cs2=[[0.0]*(n+1) for _ in range(p)]
for k in range(n):
    for j in range(p): cs1[j][k+1]=cs1[j][k]+Y[k][j]; cs2[j][k+1]=cs2[j][k]+Y[k][j]**2
def segcost(a,b):
    nn=b-a; tot=0.0
    for j in range(p):
        s1=cs1[j][b]-cs1[j][a]; s2=cs2[j][b]-cs2[j][a]; tot+=s2-s1*s1/nn
    return tot
LMIN=3
def partition(beta):                       # exact optimal partitioning (== PELT optimum)
    F=[None]*(n+1); F[0]=-beta; prev=[0]*(n+1)
    for t in range(1,n+1):
        hi=t-LMIN; cands=range(0,hi+1) if hi>=0 else [0]
        best=None; bs=0
        for a in cands:
            c=F[a]+segcost(a,t)+beta
            if best is None or c<best: best=c; bs=a
        F[t]=best; prev[t]=bs
    cps=[]; t=n
    while t>0:
        a=prev[t]
        if a>0: cps.append(a)
        t=a
    return sorted(cps)
sig2=sum(0.5*statistics.pvariance([Y[k][j]-Y[k-1][j] for k in range(1,n)]) for j in range(p))/p
beta0=p*sig2*math.log(n)                    # BIC penalty (Eq. penalty)
facs=[0.1,0.15,0.2,0.3,0.45,0.65,1.0,1.5,2.2,3.3,5.0,7.5]
curve=[(beta0*f,len(partition(beta0*f))) for f in facs]
cps=partition(beta0); K0=len(cps)
OFF={"d":7.0,"W":0.0,"A":-7.0}
for nm,arr,off in [("d",dd,OFF["d"]),("W",WW,OFF["W"]),("A",AA,OFF["A"])]:
    open(f"series_{nm}.txt","w").write("\n".join(f"{s[k]} {arr[k]+off}" for k in range(n))+"\n")
    segs=[0]+cps+[n]; lines=[]
    for a,b in zip(segs[:-1],segs[1:]):
        mu=sum(arr[a:b])/(b-a); lines.append(f"> \n{s[a]} {mu+off}\n{s[min(b,n-1)]} {mu+off}")
    open(f"step_{nm}.txt","w").write("\n".join(lines)+"\n")
open("bounds.txt","w").write("\n".join(f"> \n{s[c]} -11\n{s[c]} 11" for c in cps)+"\n")
open("labels.txt","w").write(f"{s[-1]*1.005} 7 @~d@~@-std@-\n{s[-1]*1.005} 0 W@-std@-\n{s[-1]*1.005} -7 A@-std@-\n")
open("penalty.txt","w").write("\n".join(f"{b} {K}" for b,K in curve)+"\n")
kmax=max(K for _,K in curve)+1
open("sel.txt","w").write(f"> \n{beta0} 0\n{beta0} {kmax}\n")
open("info.txt","w").write(f"{s[-1]:.0f} {beta0:.2f} {K0} {curve[0][0]:.2f} {curve[-1][0]:.2f} {kmax}\n")
PY

read SMAX BETA0 K0 BMIN BMAX KMAX < info.txt
gmt begin fig_changepoint_model pdf,png
  gmt set FONT_TITLE 11p,Helvetica FONT_ANNOT_PRIMARY 9p FONT_LABEL 10p MAP_FRAME_TYPE plain MAP_GRID_PEN_PRIMARY 0.25p,gray85
  gmt basemap -R0/$SMAX/-11/11 -JX18c/9c -BWSne+t"(b) Detected segment boundaries on the standardised series" \
      -Bxa500f100g500+l"Along-strike distance (km)" -Bya2f1+l"Standardised value (offset)"
  gmt plot bounds.txt -W0.7p,gray55,--
  gmt plot series_d.txt -W0.5p,gray65
  gmt plot series_W.txt -W0.5p,gray65
  gmt plot series_A.txt -W0.5p,gray65
  gmt plot step_d.txt -W2p,navy
  gmt plot step_W.txt -W2p,darkgreen
  gmt plot step_A.txt -W2p,darkorange
  gmt text labels.txt -F+f10p,Helvetica-Bold+jLM -N -D0.12c/0
  gmt basemap -R$BMIN/$BMAX/0/$KMAX -JX18c/4c -BWSne+t"(a) Penalty selection (number of boundaries vs @~b@~)" \
      -Bxaf+l"Penalty @~b@~" -Byaf+l"Boundaries K" -Y12.0c
  gmt plot penalty.txt -W1.5p,black
  gmt plot penalty.txt -Sc0.12c -Gblack
  gmt plot sel.txt -W1.2p,red,--
  echo "$BETA0 $K0 selected @~b@~ (BIC)" | gmt text -F+f9p,Helvetica-Bold,red+jLB -D0.15c/0.1c -N
gmt end
echo "wrote fig_changepoint_model.pdf and fig_changepoint_model.png"
