"""Direct experimental Lamb-dip traces with finite 2f detection.

The detector profile is the exact second cosine Fourier coefficient of a
Lorentz line under sinusoidal frequency modulation. Crossovers are a separate
experimental response, never transitions in an intrinsic catalog.
"""
from functools import lru_cache
import json
import numpy as np
from scipy.constants import Boltzmann,atomic_mass,c as light_speed
from .belgi import PROJECT,State
from .spin import diagonalize,dipole_kernel,dipole_kernel_q1,projection,symmetry_blocks
from .wigner import wig3jj


def lorentz_2f(x,gamma,amplitude):
    z=np.asarray(x,dtype=complex)+1j*gamma
    root=np.sqrt(z*z-amplitude*amplitude)
    root=np.where(root.imag<0,-root,root)
    # Stable branch: (root-z)/a = -a/(root+z).
    q=-amplitude/(root+z)
    raw=-2*gamma*np.imag(q*q/root)
    center=-2*gamma*(np.sqrt(gamma*gamma+amplitude*amplitude)-gamma)**2/(amplitude**2*np.sqrt(gamma*gamma+amplitude*amplitude))
    return raw/abs(center)


@lru_cache(None)
def line_geometry(ju,jl,species):
    fu=np.concatenate([np.full(u.shape[1],f) for u,f,p in symmetry_blocks(ju,species)])
    fl=np.concatenate([np.full(u.shape[1],f) for u,f,p in symmetry_blocks(jl,species)])
    pu=np.concatenate([np.full(u.shape[1],p) for u,f,p in symmetry_blocks(ju,species)])
    pl=np.concatenate([np.full(u.shape[1],p) for u,f,p in symmetry_blocks(jl,species)])
    a=np.zeros((len(fu),len(fl)));b=a.copy()
    for i,f in enumerate(fu):
        for k,g in enumerate(fl):
            if species=='E' and pu[i]==pl[k]:continue
            z=wig3jj(2*int(f),2,2*int(g),0,0,0)
            if abs(z)>1e-12:a[i,k]=1/z**2
            elif f>0:
                z=wig3jj(2*int(f),2,2*int(g),-2,2,0)
                if abs(z)>1e-12:b[i,k]=1/(f*(f+1)*z*z)
    ii,kk=np.nonzero(a+b)
    # Crossovers require a common exact hyperfine state, not just the same F.
    pair=np.array([(u,v) for u in range(len(ii)) for v in range(u+1,len(ii))
                   if ii[u]==ii[v] or kk[u]==kk[v]],int)
    return a,b,ii,kk,pair,fu,fl,pu,pl


def components(g,h):
    ju,jl,sp=g['u'].j,g['l'].j,g['species']
    eu,vu,_,_=diagonalize(ju,sp,h[g['u']]);el,vl,_,_=diagonalize(jl,sp,h[g['l']])
    a,b,i,k,pair,fu,fl,pu,pl=line_geometry(ju,jl,sp)
    d0=vu.conj().T@dipole_kernel(ju,jl,sp)@vl
    d1=vu.conj().T@dipole_kernel_q1(ju,jl,sp)@vl
    strength=(abs(d0)**2*a+abs(d1)**2*b)[i,k];strength/=strength.sum()
    freq=(eu[:,None]-el)[i,k]
    crossfreq=(freq[pair[:,0]]+freq[pair[:,1]])/2 if len(pair) else np.array([])
    crossstrength=np.sqrt(strength[pair[:,0]]*strength[pair[:,1]]) if len(pair) else np.array([])
    return freq,strength,crossfreq,crossstrength


def gaussian_2f(x,sigma,amplitude):
    theta=np.arange(64)*2*np.pi/64
    y=np.mean(np.exp(-.5*((np.asarray(x)[:,None]+amplitude*np.cos(theta))/sigma)**2)*np.cos(2*theta),axis=1)
    norm=abs(np.mean(np.exp(-.5*(amplitude*np.cos(theta)/sigma)**2)*np.cos(2*theta)))
    return y/norm


def linear_response(normal,cross,x,y,crossover=True,broad=None):
    """Exact active-set linear least squares: A>=0, 0<=C<=2 A."""
    bg=np.column_stack([np.ones(len(x)),x/max(np.ptp(x),1),
                        (x/max(np.ptp(x),1))**2 if broad is None else broad])
    candidates=[]
    if crossover and np.linalg.norm(cross)>1e-12:
        design=np.column_stack([normal,cross,bg]);v=np.linalg.lstsq(design,y,rcond=None)[0]
        if v[0]>=0 and 0<=v[1]<=2*v[0]:candidates.append((design@v,v))
    for ratio in ([0.,2.] if crossover else [0.]):
        design=np.column_stack([normal+ratio*cross,bg]);v=np.linalg.lstsq(design,y,rcond=None)[0]
        if v[0]<0:v=np.r_[0.,np.linalg.lstsq(bg,y,rcond=None)[0]]
        candidates.append((design@v,np.r_[v[0],ratio*v[0],v[1:]]))
    return min(candidates,key=lambda t:np.dot(t[0]-y,t[0]-y))


def prediction(g,sticks,nuisance,crossover=True,modulation_scale=1.):
    centre,loggamma,logmod=nuisance;gamma=np.exp(loggamma);mod=np.exp(logmod)*modulation_scale
    f,s,fc,sc=sticks;x=g['x']-centre
    normal=lorentz_2f(x[:,None]-f,gamma,mod)@s
    cross=lorentz_2f(x[:,None]-fc,gamma,mod)@sc if len(fc) else np.zeros(len(x))
    broad=None
    if g.get('thermal_background',False):
        frequency_khz=g['origin_khz'] if g['id'].startswith('Xu') else 259581350.
        sigma=frequency_khz/light_speed*np.sqrt(Boltzmann*g.get('temperature_K',293.15)/(32.026215*atomic_mass))
        broad=gaussian_2f(x,sigma,mod)
    pred,amp=linear_response(normal,cross,g['x'],g['y'],crossover,broad)
    return (pred-g['y'])*g['weight'],pred,amp


def load_profiles():
    manifest=json.loads((PROJECT/'data/processed/excited-profiles/manifest.json').read_text())
    result=[]
    for row in manifest['profiles']:
        g=dict(row);data=np.genfromtxt(PROJECT/g['file'],delimiter=',',names=True)
        g.update(x=data['offset_khz'],y=data['signal'])
        for q in ['u','l']:g[q]=State('E',g['vt'],g['j'+q],g['k'+q])
        edge=np.r_[g['y'][:max(10,len(data)//8)],g['y'][-max(10,len(data)//8):]]
        noise=float(np.median(abs(np.diff(edge)))/.9538725524)
        # A finite discrepancy floor and cap avoid pretending digitized pixels
        # are raw independent high-precision observations.
        g['signal_sigma']=max(.04,noise);g['empirical_edge_noise']=noise
        g['weight']=np.sqrt(min(1.,50/len(data)))/g['signal_sigma']
        result.append(g)
    return result
