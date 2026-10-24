"""Independent D3Q27 Taylor–Green vs. uniform-zero-flow control.

Reimplements the integer local collision LUT and streaming rules from
Taylor_green_regularity-resolved-and-solved at pinned commit 4b32817.
Run: python independent_zero_vortex_control.py
Requires: Python 3 and numpy. This is an independent NumPy implementation,
not a claim to execute the repository's original C++ program or Host8 disk.
"""
import numpy as np
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import json, csv

L=24; V=L**3; K=27; U0=.20
vel=np.array([(x,y,z) for x in (-1,0,1) for y in (-1,0,1) for z in (-1,0,1)],dtype=np.int16)
e=(vel*vel).sum(1)
pairs=defaultdict(list)
for a,b in combinations(range(K),2):
 key=tuple((vel[a]+vel[b]).tolist())+(int(e[a]+e[b]),)
 pairs[key].append((a,b))
quads=set()
for ps in pairs.values():
 for pa,pb in combinations(ps,2):
  q=tuple(sorted(pa+pb))
  if len(set(q))==4: quads.add(q)
quad=[];lut=[]
vals=np.array([[(m//3**i)%3 for i in range(4)] for m in range(81)],dtype=np.uint8)
for q in sorted(quads):
 classes=defaultdict(list)
 for code,d in enumerate(vals):
  vs=vel[list(q)]; E=e[list(q)]
  key=(int(d.sum()),)+tuple((d[:,None]*vs).sum(0).tolist())+(int(d@E),)
  classes[key].append(code)
 mapping=list(range(81))
 for key in sorted(classes):
  codes=classes[key]
  for x,y in zip(codes[::2],codes[1::2]):mapping[x]=y;mapping[y]=x
 if any(mapping[i]!=i for i in range(81)):
  quad.append(q);lut.append(mapping)
quad=np.array(quad,dtype=np.int64);lut=np.array(lut,dtype=np.uint8)
assert len(quad)==288
checks=0
for gi,q in enumerate(quad):
 for code,d in enumerate(vals):
  mapped=lut[gi,code]
  assert lut[gi,mapped]==code
  z=vals[mapped]
  assert d.sum()==z.sum()
  assert (d@vel[q]==z@vel[q]).all()
  assert d@e[q]==z@e[q]
  checks+=1
allones_code=sum(3**i for i in range(4))
changed_uniform=int(np.count_nonzero(lut[:,allones_code]!=allones_code))
print('GATES',len(quad),'INVARIANT_CHECKS',checks,'UNIFORM_1111_CHANGED_BY',changed_uniform)

ci=np.arange(V)[:,None]; shift=np.array([1,3,9,27],dtype=np.int64)
x,y,z=np.meshgrid(np.arange(L),np.arange(L),np.arange(L), indexing='ij')
kv=2*np.pi/L
ux=U0*np.sin(kv*(x+.5))*np.cos(kv*(y+.5))*np.cos(kv*(z+.5))
uy=-U0*np.cos(kv*(x+.5))*np.sin(kv*(y+.5))*np.cos(kv*(z+.5))
def round_nearest_even_momentum(u):
 q=27*u/2
 return np.clip(np.sign(q)*np.floor(np.abs(q)+0.5)*2,-4,4).astype(np.int16)
qx=round_nearest_even_momentum(ux).ravel();qy=round_nearest_even_momentum(uy).ravel()
state_tg=np.ones((V,K),dtype=np.uint8)
for moment,fp,fm,ep1,em1,ep2,em2 in [(qx,22,4,23,3,21,5),(qy,16,10,17,9,15,11)]:
 df=np.sign(moment)*(np.abs(moment)==2)
 de=np.sign(moment)*(np.abs(moment)==4)
 for plus,minus,d in [(fp,fm,df),(ep1,em1,de),(ep2,em2,de)]:
  state_tg[:,plus]=(1+d).astype(np.uint8)
  state_tg[:,minus]=(1-d).astype(np.uint8)
state_blank=np.ones((V,K),dtype=np.uint8)

def metrics(state):
 n=state.astype(np.int16)
 rho=n.sum(1)
 mom=n@vel
 twE=int(np.sum(n@e))
 U=np.zeros((V,3)); mask=rho>0;U[mask]=mom[mask]/rho[mask,None]
 u3=U.reshape(L,L,L,3)
 wx=(np.roll(u3[:,:,:,2],-1,1)-np.roll(u3[:,:,:,2],1,1))/2 - (np.roll(u3[:,:,:,1],-1,2)-np.roll(u3[:,:,:,1],1,2))/2
 wy=(np.roll(u3[:,:,:,0],-1,2)-np.roll(u3[:,:,:,0],1,2))/2 - (np.roll(u3[:,:,:,2],-1,0)-np.roll(u3[:,:,:,2],1,0))/2
 wz=(np.roll(u3[:,:,:,1],-1,0)-np.roll(u3[:,:,:,1],1,0))/2 - (np.roll(u3[:,:,:,0],-1,1)-np.roll(u3[:,:,:,0],1,1))/2
 mag2=wx**2+wy**2+wz**2
 return dict(rho_min=int(rho.min()),rho_max=int(rho.max()),net_mom=list(map(int,mom.sum(0))),twice_energy=twE,kinetic=float(np.sum(U*mom)/(2*V)),enstrophy=float(np.mean(mag2)/2),max_vorticity=float(np.sqrt(mag2.max())))

def diagnostics_init(state):
 mom=state.astype(np.int16)@vel
 actual=mom[:,0:2]/27
 target=np.stack([ux.ravel(),uy.ravel()],axis=1)
 sim=float(np.sum(actual*target)/(np.linalg.norm(actual)*np.linalg.norm(target)))
 l2=float(np.linalg.norm(actual-target)/np.linalg.norm(target))
 rms_a=float(np.sqrt(np.sum(actual**2)/V))
 rms_t=float(np.sqrt(np.sum(target**2)/V))
 return dict(similarity=sim,l2_error=l2,rms_actual=rms_a,rms_target=rms_t)

def step(state,t):
 for j in range(32):
  which=(t*32+np.arange(V)*17+j*37)%len(quad)
  q=quad[which]
  d=state[ci,q]
  codes=d.astype(np.int32)@shift
  dd=vals[lut[which,codes]]
  state[ci,q]=dd
 g=state.reshape((L,L,L,K)); nxt=np.empty_like(g)
 for k,(dx,dy,dz) in enumerate(vel):
  nxt[:,:,:,k]=np.roll(g[:,:,:,k],shift=(dx,dy,dz),axis=(0,1,2))
 return nxt.reshape(V,K)

init=metrics(state_tg); init_diag=diagnostics_init(state_tg)
print('INIT TG',json.dumps(init))
print('INIT TG FIT',json.dumps(init_diag))
print('INIT BLANK',json.dumps(metrics(state_blank)))
assert init['rho_min']==init['rho_max']==27
assert init['net_mom']==[0,0,0]
assert abs(init_diag['similarity']-0.9624369693)<1e-7
assert abs(init_diag['l2_error']-0.2793051912)<1e-7

results=[]
for label,st in [('taylor_green',state_tg.copy()),('uniform_zero_flow',state_blank.copy())]:
 for t in range(49):
  m=metrics(st);results.append(dict(control=label,tick=t,**m))
  if t in [0,1,2,3,6,12,24,48]:print(label,t,json.dumps(m))
  if t<48:st=step(st,t)
p=Path(__file__).with_name('independent_zero_vortex_results.csv')
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
print('WROTE',p)
