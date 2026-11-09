"""Uncoupled proton-spin algebra, projected onto permutation symmetry.

Four spin-1/2 protons, ordered CH3[0:3], OH[3]. Fixed total MF=0
eliminates redundant magnetic components without losing a hyperfine level.
"""
from functools import lru_cache
import itertools
import numpy as np
from .wigner_adapter import wig3jj as wigner3j, initialize

@lru_cache(None)
def spin_basis(j):
    spins=np.array(list(itertools.product([-.5,.5],repeat=4)))
    spins=spins[abs(spins.sum(axis=1))<=j]
    mj=-spins.sum(axis=1)
    lookup={tuple(s):i for i,s in enumerate(spins)}
    dim=len(spins)
    ij=np.zeros((4,dim,dim));pairs={}
    for i in range(4):
        ij[i]=np.diag(spins[:,i]*mj)
        for col,s in enumerate(spins):
            if s[i]==-.5:
                t=s.copy();t[i]=.5
                if tuple(t) in lookup:
                    row=lookup[tuple(t)]
                    v=.5*np.sqrt(j*(j+1)-mj[col]*(mj[col]-1))
                    ij[i,row,col]=ij[i,col,row]=v
    for i in range(4):
        for k in range(i+1,4):
            a=np.diag(spins[:,i]*spins[:,k])
            for col,s in enumerate(spins):
                if s[i]!=s[k]:
                    t=s.copy();t[i],t[k]=s[k],s[i]
                    a[lookup[tuple(t)],col]=.5
            pairs[i,k]=a
    cycle=np.zeros((dim,dim));exchange=np.zeros_like(cycle)
    for col,s in enumerate(spins):
        cycle[lookup[tuple(s[[2,0,1,3]])],col]=1
        exchange[lookup[tuple(s[[0,2,1,3]])],col]=1
    f2=j*(j+1)*np.eye(dim)+3*np.eye(dim)+2*sum(pairs.values())+2*ij.sum(axis=0)
    rank2={pair:1.5*(ij[pair[0]]@ij[pair[1]]+ij[pair[1]]@ij[pair[0]])-j*(j+1)*v for pair,v in pairs.items()}
    return spins,mj,ij,pairs,rank2,cycle,exchange,f2

@lru_cache(None)
def projection(j,species):
    spins,mj,ij,pairs,rank2,c,p,f2=spin_basis(j)
    d=len(spins)
    # Reuse exact saved projection basis wherever the old archive provides it.
    from pathlib import Path
    saved_path=Path(__file__).resolve().parent/'data'/'inputs/complete_frozen_projections.npz'
    if saved_path.exists():
        with np.load(saved_path,allow_pickle=False) as saved:
            key=f'{species}_{j}'
            if key+'__projection' in saved.files:
                return tuple(saved[key+suffix].copy() for suffix in ['__projection','__PI','__F2'])
    if species=='A':
        proj=(np.eye(d)+c+c@c)/3
        ev,u=np.linalg.eigh(proj)
        t=u[:,ev>.5]
        return t,t.conj().T@p@t,t.conj().T@f2@t
    # Fourier coefficients at angle alpha+2*pi*i/3 imply this convention.
    blocks=[]
    for sigma in [1,-1]:
        phase=np.exp(2j*np.pi*sigma/3)
        proj=(np.eye(d)+phase*c+phase**2*(c@c))/3
        ev,u=np.linalg.eigh((proj+proj.conj().T)/2)
        blocks.append(u[:,ev>.5])
    t=np.zeros((2*d,sum(b.shape[1] for b in blocks)),complex)
    t[:d,:blocks[0].shape[1]]=blocks[0];t[d:,blocks[0].shape[1]:]=blocks[1]
    parity=np.block([[np.zeros_like(p),p],[p,np.zeros_like(p)]])
    return t,t.conj().T@parity@t,t.conj().T@np.kron(np.eye(2),f2)@t

def project_hamiltonian(j,species,sr,ss):
    """sr[i,sigma,sigma']; ss[(i,k),sigma,sigma'] in kHz."""
    spins,mj,ij,pairs,rank2,c,p,f2=spin_basis(j)
    ns=1 if species=='A' else 2
    h=np.zeros((ns*len(spins),)*2,complex)
    for i in range(4):h+=np.kron(sr[i],ij[i])
    for pair,coeff in ss.items():h+=np.kron(coeff,rank2[pair])
    t,par,ft=projection(j,species)
    h=t.conj().T@h@t
    return (h+h.conj().T)/2

@lru_cache(None)
def symmetry_blocks(j,species):
    t,p,f2=projection(j,species)
    # F-squared and parity commute; this avoids arbitrary degenerate mixtures.
    label=f2+(.01*p if species=='E' else 0)
    vals,vec=np.linalg.eigh(label)
    blocks=[]
    for v in np.unique(np.round(vals,7)):
        use=np.flatnonzero(abs(vals-v)<1e-6)
        basis=vec[:,use]
        f=int(round((np.sqrt(1+4*np.real(basis[:,0].conj()@f2@basis[:,0]))-1)/2))
        parity=int(round(np.real(basis[:,0].conj()@p@basis[:,0]))) if species=='E' else 1
        blocks.append((basis,f,parity))
    return blocks

def diagonalize(j,species,h):
    """Resolve F and (for E) permutation-inversion parity before H."""
    rows=[];energies=[];fs=[];ps=[]
    for basis,f,parity in symmetry_blocks(j,species):
        e,u=np.linalg.eigh(basis.conj().T@h@basis)
        out=basis@u
        for ei,vi in zip(e,out.T):
            energies.append(ei);rows.append(vi);fs.append(f);ps.append(parity)
    return np.array(energies),np.array(rows).T,np.array(fs),np.array(ps)

@lru_cache(None)
def dipole_kernel(ju,jl,species,mf=0):
    """q=0 dipole in MF=0. Q-branch F=F lines need MF=1; handled below."""
    su,mu,*_=spin_basis(ju);sl,ml,*_=spin_basis(jl)
    d=np.zeros((len(su),len(sl)))
    for i,s in enumerate(su):
        for k,t in enumerate(sl):
            if np.array_equal(s,t):
                m=int(mu[i])
                d[i,k]=(-1.)**(ju-m)*wigner3j(2*ju,2,2*jl,-2*m,0,2*m)
    tu,_,_=projection(ju,species);tl,_,_=projection(jl,species)
    # Electric dipole is odd under (23)*.
    if species=='E':d=np.kron(np.diag([1.,-1.]),d)
    return tu.conj().T@d@tl

@lru_cache(None)
def dipole_kernel_q1(ju,jl,species):
    # Direct MF=1 upper eigenvectors obtained by F+ acting on MF=0.
    su,mu,*_=spin_basis(ju);sl,ml,*_=spin_basis(jl)
    allsp=np.array(list(itertools.product([-.5,.5],repeat=4)))
    spins1=allsp[abs(1-allsp.sum(axis=1))<=ju];m1=1-spins1.sum(axis=1)
    lookup={tuple(s):i for i,s in enumerate(spins1)}
    raising=np.zeros((len(spins1),len(su)))
    for c,s in enumerate(su):
        if tuple(s) in lookup:
            raising[lookup[tuple(s)],c]+=np.sqrt(ju*(ju+1)-mu[c]*(mu[c]+1))
        for i in range(4):
            if s[i]==-.5:
                t=s.copy();t[i]=.5
                if tuple(t) in lookup:raising[lookup[tuple(t)],c]+=1
    d=np.zeros((len(spins1),len(sl)))
    for i,s in enumerate(spins1):
        for k,t in enumerate(sl):
            if np.array_equal(s,t):
                d[i,k]=(-1.)**(ju-int(m1[i]))*wigner3j(2*ju,2,2*jl,-2*int(m1[i]),2,2*int(ml[k]))
    tu,_,_=projection(ju,species);tl,_,_=projection(jl,species)
    if species=='E':
        raising=np.kron(np.eye(2),raising);d=np.kron(np.diag([1.,-1.]),d)
    return (raising@tu).conj().T@d@tl

def transition_components(ju,jl,species,hu,hl,return_indices=False):
    """Relative unpolarized strengths from q=0/MF=0 and q=1/MF'=1.

    The second route covers F'=F'', invisible in the first route. Both
    use reduced F matrix elements and the same spin-independent dipole.
    """
    eu,vu,fu,pu=diagonalize(ju,species,hu)
    el,vl,fl,pl=diagonalize(jl,species,hl)
    d0=vu.conj().T@dipole_kernel(ju,jl,species)@vl
    d1=vu.conj().T@dipole_kernel_q1(ju,jl,species)@vl
    result=[]
    for i in range(len(eu)):
        for k in range(len(el)):
            a=wigner3j(2*int(fu[i]),2,2*int(fl[k]),0,0,0)
            if abs(a)>1e-12:strength=abs(d0[i,k]/a)**2
            elif fu[i]>0:
                a=wigner3j(2*int(fu[i]),2,2*int(fl[k]),-2,2,0)
                strength=abs(d1[i,k]/(np.sqrt(fu[i]*(fu[i]+1))*a))**2 if abs(a)>1e-12 else 0
            else:strength=0
            if strength>1e-12:
                if species=='E' and pu[i]==pl[k]:
                    raise ArithmeticError('Electric-dipole component violates total G6 A1 <-> A2 selection')
                row=(eu[i]-el[k],strength,int(fu[i]),int(fl[k]),int(pu[i]),int(pl[k]))
                result.append(row+(i,k) if return_indices else row)
    return np.array(result)
