"""Rebuild L26 coefficient operators from corrected theory tensors and BELGI waves."""
import json,ast
from pathlib import Path
import numpy as np
from methanol_hyperfine.tensors import Contractions,fourier_tensor,RHO
from methanol_hyperfine.spin import project_hamiltonian,projection
from methanol_hyperfine.belgi import export_wavefunctions
from .forward import DATA,load_problem
THEORY=json.loads((DATA/'theory_tensors.json').read_text())
PUBLISHED=dict(spin_rotation=THEORY['M'],spin_torsion_total=THEORY['total_d'])
ST_TERMS=[('CH3','a',0),('CH3','b',0),('OH','a',0),('OH','b',0),('CH3','a',1)]
EXTRA_ST=[]
SR_TERMS=[(g,c) for g,t in PUBLISHED['spin_rotation'].items() for c in t]
MASKS={key:lambda n,g:n<=2 for key in ['SR','ST','DD']}
def frame_parameters(frame='paper-vector'):
    if frame!='paper-vector':raise ValueError('L26 uses the paper-vector frame')
    physical=json.loads((DATA.parents[2]/'methanol_zeeman/data/inputs/physical_parameters.json').read_text())['physical_tensors']
    lam=np.array(physical['lambda_RAM']);return lam,1/(1-RHO*lam[0])
def dipolar_fourier(order=12):
    result={}
    for name,item in THEORY['DD'].items():
        ft={}
        for n,(co,si) in enumerate(zip(item['cosine_coefficients_khz'],item['sine_coefficients_khz'])):
            ft.update(fourier_tensor(co,si,n))
        result[ast.literal_eval(name)]=ft
    return result,{}
def prepare_masked_state(w,frame='paper-vector',extended=False,dd_order=12):
    terms=ST_TERMS+(EXTRA_ST if extended else []);n_st=len(terms)
    j=w.state.j;species=w.state.species;ns=1 if species=='A' else 2
    dim=projection(j,species)[0].shape[1]
    fixed=np.zeros((dim,dim),complex);basis=np.zeros((n_st+20,dim,dim),complex)
    if not j:return fixed,basis
    c=Contractions(w);jj=j*(j+1);lam,factor=frame_parameters(frame)
    blank=np.zeros((4,ns,ns),complex)
    sr_fixed=blank.copy();sr_basis=np.zeros((n_st+18,4,ns,ns),complex)
    for proton in range(4):
        group='CH3' if proton<3 else 'OH';phase=2*np.pi*proton/3 if proton<3 else 0.
        # Total d, keeping previously calculated/fitted non-selected harmonics explicit.
        for axis,a in enumerate('abc'):
            for harmonic,value in enumerate(PUBLISHED['spin_torsion_total'][group][a]):
                if not MASKS.get('ST',lambda n,g:True)(harmonic,group):continue
                term=(group,a,harmonic)
                op=factor*c.spin_torsion(axis,harmonic,phase,axis==2)/jj
                if term in terms:sr_basis[terms.index(term),proton]+=op
                else:sr_fixed[proton]+=value*op
        # -M plus the -M.lambda contribution to d' must vary together.
        for component,values in PUBLISHED['spin_rotation'][group].items():
            a,b=map('abc'.index,component);operator=np.zeros((ns,ns),complex)
            for harmonic,value in enumerate(values):
                if not value or not MASKS.get('SR',lambda n,g:True)(harmonic,group):continue
                cosine=np.zeros((3,3));sine=np.zeros((3,3))
                (sine if (a==2)!=(b==2) else cosine)[a,b]=value
                operator-=c.tensor(fourier_tensor(cosine,sine,harmonic,phase))/jj
                if lam[b]:operator-=factor*lam[b]*value*c.spin_torsion(a,harmonic,phase,a==2)/jj
            sr_fixed[proton]+=operator
            sr_basis[n_st+SR_TERMS.index((group,component)),proton]+=operator
    fixed=project_hamiltonian(j,species,sr_fixed,{})
    for k,sr in enumerate(sr_basis):basis[k]=project_hamiltonian(j,species,sr,{})
    dd,_=dipolar_fourier(dd_order)
    dd={pair:{n:v for n,v in coeff.items() if MASKS.get('DD',lambda n,g:True)(abs(n),'CH3' if pair[1]<3 else 'OH')} for pair,coeff in dd.items()}
    coefficients={pair:2*c.tensor(coeff)/(jj*(2*j-1)*(2*j+3)) for pair,coeff in dd.items()}
    for k,oh in enumerate([False,True]):
        op=project_hamiltonian(j,species,blank,{pair:v for pair,v in coefficients.items() if (pair[1]==3)==oh})
        fixed+=op;basis[n_st+18+k]=op
    return fixed,basis

def verify_operators():
    p=load_problem(wavefunctions="reference");d=p.d;waves=export_wavefunctions(d.prepared)
    errors={}
    for state,w in waves.items():
        fixed,basis=prepare_masked_state(w,d.frame);basis=basis*d.coordinate_to_fraction[:,None,None]
        reference=d.prepared[state]
        errors[state.key]=max(float(np.max(abs(fixed-reference[0]))),float(np.max(abs(basis-reference[1]))))
    maximum=max(errors.values())
    if maximum>=1e-10:raise ArithmeticError('Operator regeneration mismatch: '+str(maximum))
    return dict(passed=True,states=len(errors),max_operator_error_khz=maximum)
if __name__=='__main__':print(json.dumps(verify_operators(),indent=2))
