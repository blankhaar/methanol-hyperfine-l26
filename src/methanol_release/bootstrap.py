"""Conditional local-quadratic uncertainty draws in the L26 seven-mode basis.

These ranges condition on fixed tensors, assignments, weighting, basis and
assumed noise. They are not calibrated experimental confidence intervals.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,lsq_linear
from .forward import DATA as FIT_DATA, load_problem
from methanol_zeeman.output_paths import validate_output_directory
DATA=Path(__file__).parent/'data/uncertainty'
SELECTIVE_ROWS=6432
SEED=20261019
class QuadraticBootstrap:
    def __init__(self, wavefunctions="native"):
        if wavefunctions not in ('native','reference'):
            raise ValueError('Choose native or reference wavefunctions')
        final=np.load(FIT_DATA/'endpoint.npz',allow_pickle=False)
        basis=np.load(FIT_DATA/'basis.npz',allow_pickle=False)
        values=load_problem().evaluate(final['x'],True) if wavefunctions=='native' else dict(jacobian=final['jacobian'],residual=final['residual'])
        self.wavefunction_source=wavefunctions
        self.M_basis=basis['M_basis'];raw=values['jacobian']
        self.J=np.column_stack([raw[:,:4],raw[:,4:22]@self.M_basis,raw[:,22:]])
        self.r=values['residual'];self.ra=self.r[:SELECTIVE_ROWS]
        self.s=1/np.maximum(np.linalg.norm(self.J,axis=0),1e-10)
        self.Js=self.J[:SELECTIVE_ROWS]*self.s;self.Jl=self.J[SELECTIVE_ROWS:]*self.s
        self.lo=np.r_[final['lower'][:4]-final['x'][:4],np.full(7,-np.inf),final['lower'][22:]-final['x'][22:]]/self.s
        self.hi=np.r_[final['upper'][:4]-final['x'][:4],np.full(7,np.inf),final['upper'][22:]-final['x'][22:]]/self.s
        keep=np.linalg.norm(self.Js,axis=0)>1e-10
        self.A=self.Js[:,keep];self.alo=self.lo[keep];self.ahi=self.hi[keep]
        anchor=lsq_linear(self.A,-self.ra,bounds=(self.alo,self.ahi),method='bvls',tol=1e-12,max_iter=500)
        ceiling=json.loads((FIT_DATA/'results.json').read_text())['selective_ceiling']
        self.allowance=ceiling-anchor.cost
        self.Ha=self.A.T@self.A;self.Hs=self.Js.T@self.Js;self.Hl=self.Jl.T@self.Jl
        archive=np.load(DATA/'draws.npz',allow_pickle=False)
        self.D=archive['derivatives']
        self.centers=archive['centers'];self.labels=archive['labels']
        if wavefunctions=='reference' and abs(self.allowance-4.640706911811662)>1e-9:raise ArithmeticError('Selective anchor mismatch')
    def solve(self, noise: np.ndarray) -> tuple[dict, np.ndarray | None]:
            """Solve both bounded SLSQP stages for one residual perturbation.

            The SAME selective noise enters the anchor and reduced stages. The
            final guard is the fitted perturbed anchor Q plus the fixed allowance.
            Failure to find a feasible point is not, by itself, proof of infeasibility.
            """
            e = np.asarray(noise, dtype=float)
            if e.shape != self.r.shape or not np.isfinite(e).all():
                raise ValueError("Noise must be a finite length-10142 residual vector")
            started = time.monotonic()
            ys, yl = self.r[:SELECTIVE_ROWS] - e[:SELECTIVE_ROWS], self.r[SELECTIVE_ROWS:] - e[SELECTIVE_ROWS:]
            ya = self.ra - e[:SELECTIVE_ROWS]
            ga, gs, gl = self.A.T @ ya, self.Js.T @ ys, self.Jl.T @ yl
            qa, qs, ql = .5 * ya @ ya, .5 * ys @ ys, .5 * yl @ yl

            def anchor_objective(z):
                return .5 * z @ self.Ha @ z + ga @ z, self.Ha @ z + ga

            anchor = minimize(anchor_objective, np.zeros(self.A.shape[1]), jac=True,
                              bounds=list(zip(self.alo, self.ahi)), method="SLSQP",
                              options={"maxiter": 1000, "ftol": 1e-10})
            anchor_q = qa + anchor.fun
            ceiling = anchor_q + self.allowance

            def objective(z):
                return (.5 * z @ self.Hl @ z + gl @ z) / 100., (self.Hl @ z + gl) / 100.

            def constraint(z):
                return (ceiling - qs - gs @ z - .5 * z @ self.Hs @ z) / 10.

            def constraint_jacobian(z):
                return -(self.Hs @ z + gs) / 10.

            initial = np.zeros(self.J.shape[1])
            if constraint(initial) < 0:
                def selective_objective(z):
                    return .5 * z @ self.Hs @ z + gs @ z, self.Hs @ z + gs

                feasibility = minimize(selective_objective, initial, jac=True,
                                       bounds=list(zip(self.lo, self.hi)), method="SLSQP",
                                       options={"maxiter": 1000, "ftol": 1e-10})
                initial = feasibility.x
                if constraint(initial) < -1e-6:
                    return {
                        "valid": False, "status": "feasibility_search_failed",
                        "reason": "Original SLSQP feasibility search did not reach the perturbed anchor ceiling",
                        "slack": float(constraint(initial) * 10),
                        "anchor_success": bool(anchor.success),
                        "feasibility_success": bool(feasibility.success),
                        "seconds": time.monotonic() - started,
                    }, None

            fitted = minimize(objective, initial, jac=True, bounds=list(zip(self.lo, self.hi)),
                              constraints=[{"type": "ineq", "fun": constraint, "jac": constraint_jacobian}],
                              method="SLSQP", options={"maxiter": 1500, "ftol": 1e-10})
            slack = constraint(fitted.x) * 10
            bound_error = max(float(np.max(self.lo - fitted.x)), float(np.max(fitted.x - self.hi)), 0.)
            delta = self.s * fitted.x
            valid = bool(anchor.success and fitted.success and slack >= -1e-5 and bound_error < 1e-7
                         and np.isfinite(delta).all())
            result = dict(
                valid=valid, status="valid" if valid else "solver_failure",
                anchor_success=bool(anchor.success), stage_success=bool(fitted.success),
                anchor_message=str(anchor.message), stage_message=str(fitted.message),
                anchor_Q=float(anchor_q),
                selective_Q=float(qs + gs @ fitted.x + .5 * fitted.x @ self.Hs @ fitted.x),
                global_Q=float(ql + gl @ fitted.x + .5 * fitted.x @ self.Hl @ fitted.x),
                slack=float(slack), bound_violation=bound_error,
                active_response_count=int(np.sum((abs(fitted.x[13:] - self.lo[13:]) < 1e-5) |
                                                 (abs(fitted.x[13:] - self.hi[13:]) < 1e-5))),
                seconds=time.monotonic() - started,
            )
            return result, delta

def noise_stream():
    rng=np.random.default_rng(SEED)
    for sigma,n in ((1.,512),(.5,128),(2.,128)):
        for _ in range(n):
            e=rng.normal(size=10142);e[:99]*=sigma
            yield sigma,e

def replay(indices=(0,512,640)):
    model=QuadraticBootstrap(wavefunctions="reference");reference=np.load(DATA/'draws.npz',allow_pickle=False)
    status,zero=model.solve(np.zeros(10142))
    if not status['valid']:raise ArithmeticError('Zero-noise solve failed')
    checks=[]
    for i,(sigma,noise) in enumerate(noise_stream()):
        if i not in indices:continue
        status,delta=model.solve(noise)
        if not status['valid']:raise ArithmeticError('Draw failed: '+str(i))
        error=float(np.max(abs(model.D@(delta-zero)[:13]-reference['element_displacements'][i])))
        if error>=1e-8:raise ArithmeticError('Draw replay mismatch: '+str(error))
        checks.append(dict(index=i,table_sigma=sigma,max_coefficient_error=error))
    return dict(passed=True,fit='L26',linearization='archived conditional-draw reference',same_basis_anchor_allowance=model.allowance,checks=checks,archived_draws=768,
        uncertainty='conditional local-quadratic ranges; not calibrated confidence intervals')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--replay',action='store_true');p.add_argument('--run',action='store_true');p.add_argument('--output',type=Path)
    a=p.parse_args(argv)
    if not a.run:print(json.dumps(replay(),indent=2));return 0
    if a.output is None:p.error('--run requires --output')
    out=validate_output_directory(a.output,('draws.npz','diagnostics.json'));out.mkdir(parents=True,exist_ok=True)
    m=QuadraticBootstrap();status,zero=m.solve(np.zeros(10142))
    if not status['valid']:raise ArithmeticError('Zero-noise solve failed')
    values=[];coords=[];sigmas=[];diagnostics=[]
    for sigma,e in noise_stream():
        status,d=m.solve(e);diagnostics.append(status);sigmas.append(sigma)
        if not status['valid']:raise ArithmeticError('Invalid draw; no intervals emitted')
        coords.append(d-zero);values.append(m.D@(d-zero)[:13])
    np.savez_compressed(out/'draws.npz',element_displacements=values,coordinate_displacements=coords,
        table_sigma=sigmas,valid=np.ones(768,bool),centers=m.centers,labels=m.labels,derivatives=m.D,
        zero_noise_coordinate_displacement=zero)
    (out/'diagnostics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    print(json.dumps(dict(passed=True,draws=768,output=str(out))));return 0
if __name__=='__main__':raise SystemExit(main())
