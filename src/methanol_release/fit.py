"""Replay and optimize L26 in its declared affine seven-mode M subspace."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from .forward import DATA,load_problem,reduced_transform
from methanol_zeeman.output_paths import validate_output_directory

def audit(derivatives=False, *, wavefunctions="native"):
    a=np.load(DATA/'endpoint.npz',allow_pickle=False)
    b=np.load(DATA/'basis.npz',allow_pickle=False)
    endpoint=json.loads((DATA/'results.json').read_text())
    p=load_problem(wavefunctions=wavefunctions);c=p.evaluate(a['x'],derivatives);r=c['residual']
    e=float(np.max(abs(r-a['residual'])))
    mode_error=float(np.max(abs(b['M_offset']+b['M_basis']@endpoint['M_amplitudes']-a['x'][4:22])))
    result=dict(fit='L26',wavefunction_provider=wavefunctions,forward_residual_recomputed=True,residual_count=len(r),stored_coordinates=len(a['x']),independent_coordinates=127,
        max_residual_error=e,affine_M_error_khz=mode_error,
        Q_total=float(.5*r@r),Q_selective=float(.5*r[:6432]@r[:6432]),Q_low=float(.5*r[6432:]@r[6432:]),
        selective_ceiling=endpoint['selective_ceiling'])
    if derivatives:
        result['max_jacobian_error']=float(np.max(abs(c['jacobian']-a['jacobian'])))
        result['max_scaled_jacobian_error']=float(np.max(abs(c['jacobian']-a['jacobian'])/np.maximum(1.,abs(a['jacobian']))))
    derivative_failed = (result.get('max_jacobian_error',0)>=1e-10 if wavefunctions=='reference'
                         else result.get('max_scaled_jacobian_error',0)>=1e-7)
    result['jacobian_comparison'] = ('strict archived-operator replay, absolute 1e-10' if wavefunctions=='reference'
        else 'native operator reconstruction: 1e-7 relative to max(1,abs(reference)); finite-difference roundoff reported separately')
    if e>=1e-10 or mode_error>=1e-10 or derivative_failed:
        raise ArithmeticError('L26 replay failed: '+str(result))
    if result['Q_selective']>endpoint['selective_ceiling']+1e-5:raise ArithmeticError('Selective constraint violated')
    result['passed']=True
    return result

def optimize(output,max_iterations=1200):
    output=validate_output_directory(output,('results.json','optimizer.npz'))
    output.mkdir(parents=True,exist_ok=True)
    a=np.load(DATA/'endpoint.npz',allow_pickle=False);endpoint=json.loads((DATA/'results.json').read_text())
    p=load_problem();x0=a['x'];T=reduced_transform();J=p.jacobian(x0)@T
    scale=np.minimum(100.,1/np.maximum(np.linalg.norm(J,axis=0),.01));T=T*scale
    lo=np.r_[a['lower'][:4]-x0[:4],np.full(7,-np.inf),a['lower'][22:]-x0[22:]]/scale
    hi=np.r_[a['upper'][:4]-x0[:4],np.full(7,np.inf),a['upper'][22:]-x0[22:]]/scale
    limit=endpoint['selective_ceiling'];oscale=max(1,endpoint['global_profile_cost']/20)
    def objective(y):
        c=p.evaluate(x0+T@y,True);r=c['residual'][6432:]
        return .5*r@r/oscale,T.T@c['jacobian'][6432:].T@r/oscale
    def guard(y):
        c=p.evaluate(x0+T@y,True);r=c['residual'][:6432]
        return (limit-.5*r@r),-T.T@c['jacobian'][:6432].T@r
    fitted=minimize(objective,np.zeros(127),jac=True,method='SLSQP',bounds=list(zip(lo,hi)),
        constraints=[dict(type='ineq',fun=lambda y:guard(y)[0],jac=lambda y:guard(y)[1])],
        options=dict(maxiter=max_iterations,ftol=1e-8))
    x=x0+T@fitted.x;c=p.evaluate(x,True);r=c['residual']
    report=dict(fit='L26',wavefunction_provider='native BELGI',success=bool(fitted.success and guard(fitted.x)[0]>=-1e-5),optimizer_success=bool(fitted.success),
        message=str(fitted.message),iterations=int(fitted.nit),evaluations=int(fitted.nfev),
        Q_total=float(.5*r@r),Q_selective=float(.5*r[:6432]@r[:6432]),Q_low=float(.5*r[6432:]@r[6432:]),
        selective_ceiling=limit,coordinate_values=x.tolist(),parameters=p.unpack(x).tolist(),
        start='packaged L26',M_basis='packaged affine seven-mode basis',replaces_packaged_fit=False)
    (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(output/'optimizer.npz',x=x,residual=r,jacobian=c['jacobian'],lower=a['lower'],upper=a['upper'])
    return report

def main(argv=None):
    q=argparse.ArgumentParser(description=__doc__);q.add_argument('--replay',action='store_true');q.add_argument('--derivatives',action='store_true');q.add_argument('--run',action='store_true');q.add_argument('--output',type=Path);q.add_argument('--max-iterations',type=int,default=1200);q.add_argument('--reference',action='store_true',help='Explicit exact archived-operator replay; not the default native calculation')
    a=q.parse_args(argv)
    if a.run:
        if a.reference:q.error('--reference is only for archived replay')
        if not a.output:q.error('--run requires a new --output directory')
        result=optimize(a.output,a.max_iterations)
    else:result=audit(a.derivatives,wavefunctions='reference' if a.reference else 'native')
    print(json.dumps(result,indent=2));return 0 if result.get('passed',result.get('success',False)) else 2
if __name__=='__main__':sys.exit(main())
