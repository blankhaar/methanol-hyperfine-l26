"""Spatial contractions in the BELGI basis; (a,b,c)=(z,x,y), exp(+i n alpha)."""
import numpy as np
RHO=.8102062229580776
class Contractions:
    def __init__(self,w):
        self.w=w;self.j=j=w.state.j;self.jj=j*(j+1)
        k=np.arange(-j,j+1);t=w.ktronc;self.ms=np.arange(-3*t-1,3*t+2)
        ns=1 if w.state.species=='A' else 2
        self.phi=np.zeros((ns,2*j+1,len(self.ms)),complex)
        sigma=int(w.state.species=='E')
        for n in range(-t,t+1):self.phi[0,:,3*n+sigma-self.ms[0]]=w.coefficients[:,n+t]
        if ns==2:self.phi[1]=self.phi[0,::-1,::-1]*((-1.)**(j-k))[:,None]
        jp=np.zeros((2*j+1,2*j+1),complex)
        for col,kk in enumerate(k):
            if col>0:jp[col-1,col]=np.sqrt(jj_value(j)-kk*(kk-1))
        self.J=[np.diag(k),.5*(jp+jp.T),(.5/1j)*(jp-jp.T)]
        self.cache={}

    def quadratic(self,a,b,n):
        key=a,b,n
        if key in self.cache:return self.cache[key]
        mlen=len(self.ms)
        if n>=0:left=self.phi[:,:,n:];right=self.phi[:,:,:mlen-n]
        else:left=self.phi[:,:,:mlen+n];right=self.phi[:,:,-n:]
        r=self.J[a]@self.J[b]
        rr,cc=np.nonzero(r)
        value=np.einsum('ark,r,brk->ab',left[:,rr,:].conj(),r[rr,cc],right[:,cc,:],optimize=True)
        self.cache[key]=value
        return value

    def tensor(self,coeff):
        result=np.zeros((self.phi.shape[0],)*2,complex)
        for n,t in coeff.items():
            for a,b in zip(*np.nonzero(abs(t)>1e-10)):
                result+=t[a,b]*self.quadratic(a,b,n)
        return (result+result.conj().T)/2

    def momentum(self,axis,n):
        """Hermitian symmetrization of J_axis exp(in alpha) p_alpha."""
        key=('p',axis,n)
        if key in self.cache:return self.cache[key]
        length=len(self.ms)
        if n>=0:
            left=self.phi[:,:,n:];right=self.phi[:,:,:length-n]
            m=(self.ms[n:]+self.ms[:length-n])/2
        else:
            left=self.phi[:,:,:length+n];right=self.phi[:,:,-n:]
            m=(self.ms[:length+n]+self.ms[-n:])/2
        rr,cc=np.nonzero(self.J[axis])
        value=np.einsum('ark,r,k,brk->ab',left[:,rr,:].conj(),self.J[axis][rr,cc],m,right[:,cc,:],optimize=True)
        self.cache[key]=value
        return value

    def spin_torsion(self,axis,harmonic,phase=0.,sine=False):
        if harmonic==0:coeff={0:1.}
        else:
            z=(.5/1j if sine else .5)*np.exp(1j*harmonic*phase)
            coeff={harmonic:z,-harmonic:np.conj(z)}
        result=np.zeros((self.phi.shape[0],)*2,complex)
        for n,z in coeff.items():
            result+=z*(self.momentum(axis,n)+RHO*.5*(self.quadratic(axis,0,n)+self.quadratic(0,axis,n)))
        return (result+result.conj().T)/2

def jj_value(j):return j*(j+1)

def fourier_tensor(cosine,sine,harmonic,phase=0.):
    if harmonic==0:return {0:np.asarray(cosine,dtype=complex)}
    value=(np.asarray(cosine)-1j*np.asarray(sine))/2*np.exp(1j*harmonic*phase)
    return {harmonic:value,-harmonic:value.conj()}

def state_hamiltonian(prepared,parameters):
    fixed,basis=prepared
    return fixed+np.einsum('i,ijk->jk',parameters,basis)
