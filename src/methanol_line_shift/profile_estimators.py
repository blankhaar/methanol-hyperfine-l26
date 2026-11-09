"""Noise-free line estimators, with closed-form full-real-axis inner products.

All frequency inputs and returned positions are in the same units. Width W is
component ABSORPTION FWHM. Uniform L2, known zero baseline, same response kernel,
free positive amplitude/center/absorption width. Derivative orders label ideal
zero-dither linear response. Gaussian: 0f absorption; Lorentz: 0f absorption,
1f dispersion derivative (same autocorrelation as absorption derivative), and
2f absorption derivative. No population, maser, or optical excitation model.
"""
import math
import numpy as np
from scipy.optimize import root, minimize, brentq, minimize_scalar
from scipy.special import eval_hermitenorm

SHAPES=('Gaussian','Lorentz','dispersion_1f','absorption_2f')

def overlap_and_derivatives(d,t,base,shape):
    # Source widths base and fit widths t in sigma or gamma, d=f_i-center.
    if shape=='Gaussian':
        v=base*base+t*t
        c=np.exp(-d*d/(2*v))/np.sqrt(2*np.pi*v)
        dc=c*d/v
        dt=c*t*(d*d/v**2-1/v)
        norm=1/(2*np.sqrt(np.pi)*t); exponent=1
    else:
        order={'Lorentz':0,'dispersion_1f':1,'absorption_2f':2}[shape]
        m=2*order+1;z=base+t-1j*d; fac=math.factorial(2*order)/np.pi
        c=fac*np.real(z**(-m))
        dc=fac*np.real(-m*1j*z**(-m-1))
        dt=fac*np.real(-m*z**(-m-1))
        norm=fac/(2*t)**m;exponent=m
    return c,dc,dt,norm,exponent

def fitted_center(offset,strength,width,shape='Gaussian'):
    f=np.asarray(offset,dtype=float);s=np.asarray(strength,dtype=float);s=s/s.sum()
    # Dimensionless calculation protects sub-Hz pulls at broad widths.
    centre0=float(s@f);d=(f-centre0)/width
    base=1/np.sqrt(8*np.log(2)) if shape=='Gaussian' else .5
    def stats(z):
        c,dc,dt,norm,m=overlap_and_derivatives(d-z[0],np.exp(z[1]),base,shape)
        return s@c,s@dc,s@dt,norm,m
    def equation(z):
        with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
            C,Cc,Ct,N,m=stats(z);t=np.exp(z[1])
            value=np.array([Cc/C,2*t*Ct/C+m])
        return np.nan_to_num(value,nan=1e8,posinf=1e8,neginf=-1e8)
    def objective(z):
        C,Cc,Ct,N,m=stats(z)
        return -C*C/N if C>0 else 1e6
    guesses=[0.,float(s@d),float(d[np.argmax(s)]),float(d.min()),float(d.max())]
    candidates=[]
    for c0 in sorted(set(guesses)):
        sol=root(equation,[c0,np.log(base)],tol=1e-11)
        z=sol.x
        if np.all(np.isfinite(z)) and abs(z[0])<max(abs(d))+3 and abs(z[1])<8 and np.linalg.norm(equation(z))<1e-8 and stats(z)[0]>0:
            candidates.append(z)
    if not candidates:raise RuntimeError((shape,width,'No stationary fit'))
    z=min(candidates,key=objective)
    # The tested widths are unresolved. Verify local maximum of squared projection.
    eps=1e-5
    assert all(objective(z)<=objective(z+np.eye(2)[k]*eps)+1e-11 for k in range(2))
    C,_,_,N,_=stats(z)
    return dict(center=centre0+z[0]*width,fitted_absorption_FWHM=np.exp(z[1])/base*width,amplitude=C/N,stationarity_norm=float(np.linalg.norm(equation(z))))

def kernel(x,width,shape):
    if shape=='Gaussian':
        sigma=width/np.sqrt(8*np.log(2));return np.exp(-.5*(x/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
    gamma=width/2
    if shape=='Lorentz':return gamma/(np.pi*(x*x+gamma*gamma))
    if shape=='dispersion_1f':return (gamma*gamma-x*x)/(np.pi*(x*x+gamma*gamma)**2)
    if shape=='absorption_2f':return 2*gamma*(3*x*x-gamma*gamma)/(np.pi*(x*x+gamma*gamma)**3)
    raise ValueError(shape)

def absorption_fwhm(offset,strength,width,shape='Gaussian'):
    """FWHM of the composite ABSORPTION profile, not the derivative trace."""
    shape='Gaussian' if shape=='Gaussian' else 'Lorentz'
    f=np.array(offset);s=np.array(strength);s=s/s.sum()
    def y(x):return float(s@kernel(x-f,width,shape))
    peak=minimize_scalar(lambda x:-y(x),bounds=(float(f.min()-width),float(f.max()+width)),method='bounded',options={'xatol':1e-12}).x
    half=y(peak)/2
    left=brentq(lambda x:y(x)-half,float(f.min()-10*width),peak,xtol=1e-12)
    right=brentq(lambda x:y(x)-half,peak,float(f.max()+10*width),xtol=1e-12)
    return right-left

def component_width_for_composite(offset,strength,composite_width,shape='Gaussian'):
    return brentq(lambda w:absorption_fwhm(offset,strength,w,shape)-composite_width,composite_width*.5,composite_width,xtol=1e-10)

def broad_limit(offset,strength,width,shape='Gaussian'):
    f=np.array(offset);s=np.array(strength);s=s/s.sum();mu3=float(s@(f-s@f)**3)
    coef={'Gaussian':2*np.log(2),'Lorentz':2.,'dispersion_1f':5.,'absorption_2f':28/3}[shape]
    return -coef*mu3/width**2

def component_width_for_fitted_fwhm(offset,strength,fitted_width,shape='Gaussian'):
    """Choose component width to match a published single-profile fit width."""
    return brentq(lambda w:fitted_center(offset,strength,w,shape)['fitted_absorption_FWHM']-fitted_width,fitted_width*.5,fitted_width,xtol=1e-10)
