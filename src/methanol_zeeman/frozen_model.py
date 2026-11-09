"""Forward contraction of L26 physical tensors. No optimization.

The supplied dprime and M are already in native conventions. Only magnetic G
gets the separately adopted single-c conversion in catalog.contract_magnetic.
"""
import json,ast
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from .spin import project_hamiltonian,projection
from .contractions import Contractions,fourier_tensor,RHO
O=Path(__file__).resolve().parent/'data'
parameters=json.loads((O/'inputs/physical_parameters.json').read_text())
physical=parameters['physical_tensors'];factor=1/(1-RHO*physical['lambda_RAM'][0])
DD=json.loads((O/'inputs/DD_fourier_tensors.json').read_text())
def wave(key,coeff):
 """Minimal immutable-wavefunction view expected by spatial contractions."""
 x=key.split('_');return SimpleNamespace(state=SimpleNamespace(j=int(x[2]),species=x[0]),ktronc=10,coefficients=coeff)
def hyperfine(key,coeff):
 """Contract frozen d-prime, M and DD tensors to the physical spin subspace.

 Returns H/h in kHz. This function has no optimizer, no free parameters and no
 inter-level mixing. The d-prime series already include the M-lambda term;
 nuclear subtraction is a reporting operation and must not be applied here.
 """
 w=wave(key,coeff);j=w.state.j;ns=1 if w.state.species=='A' else 2
 if not j: return np.zeros((projection(j,w.state.species)[0].shape[1],)*2)
 con=Contractions(w);jj=con.jj;sr=np.zeros((4,ns,ns),complex)
 for proton in range(4):
  group='CH3' if proton<3 else 'OH';phase=2*np.pi*proton/3 if proton<3 else 0.
  for a,axis in enumerate('abc'):
   for n,value in enumerate(physical['effective_d_prime'][group][axis]):
    if value:sr[proton]+=value*factor*con.spin_torsion(a,n,phase,a==2)/jj
  for pair,series in physical['M'][group].items():
   a,b=map('abc'.index,pair)
   for n,value in enumerate(series):
    if not value:continue
    cosine=np.zeros((3,3));sine=cosine.copy()
    (sine if (a==2)!=(b==2) else cosine)[a,b]=value
    sr[proton]-=con.tensor(fourier_tensor(cosine,sine,n,phase))/jj
 ss={}
 for pair,item in DD.items():
  ft={}
  for n,(cosine,sine) in enumerate(zip(item['cosine_coefficients_khz'],item['sine_coefficients_khz'])):
   ft.update(fourier_tensor(cosine,sine,n))
  ss[ast.literal_eval(pair)]=2*con.tensor(ft)/(jj*(2*j-1)*(2*j+3))
 return project_hamiltonian(j,w.state.species,sr,ss)
if __name__=='__main__':
 from methanol_hyperfine.belgi import State, export_wavefunctions
 ref=np.load(O/'inputs/reference_hamiltonians.npz');errors={}
 keys=np.load(O/'inputs/available_wavefunctions.npz').files
 waves={s.key:w.coefficients for s,w in export_wavefunctions([State.from_key(k) for k in keys]).items()}
 for key in ref.files:
  if key in waves:
   errors[key]=float(np.max(abs(hyperfine(key,waves[key])-ref[key])))
 print('Compared states',len(errors),'maximum Hamiltonian error/kHz',max(errors.values()))
 print(sorted(errors.items(),key=lambda x:-x[1])[:10])
 # Standalone inspection is read-only. Use the public CLI for explicit outputs.
 assert max(errors.values())<1e-9
