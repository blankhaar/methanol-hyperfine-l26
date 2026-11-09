"""Independent native/reference and derivative checks; outputs outside source."""
import argparse, hashlib, json, os, sys, time
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--mode',choices=['waves','derivatives','bootstrap'],required=True)
a=p.parse_args();sys.path.insert(0,str(a.source/'src'))
from methanol_hyperfine.belgi import State,export_wavefunctions,PROJECT
from methanol_hyperfine.native_belgi import configuration
start=time.monotonic();report={'mode':a.mode,'configuration':configuration()[1]}

if a.mode=='waves':
    results=[]
    for k in [10,12,14]:
        meta=json.loads((PROJECT/f'data/wavefunctions/ktronc-{k}.json').read_text())
        states=[State.from_key(s) for s in meta['energies']]
        fresh=export_wavefunctions(states,ktronc=k,force=True)
        reference=export_wavefunctions(states,ktronc=k,backend='reference')
        energy=max(abs(fresh[s].energy-reference[s].energy) for s in states)
        coef=max(min(np.max(abs(fresh[s].coefficients-reference[s].coefficients)),
                     np.max(abs(fresh[s].coefficients+reference[s].coefficients))) for s in states)
        norm=max(abs(fresh[s].norm-reference[s].norm) for s in states)
        item=dict(ktronc=k,states=len(states),max_energy_error_cm_inverse=float(energy),
            max_coefficient_error_up_to_global_sign=float(coef),max_norm_error=float(norm),
            passed=bool(energy<1e-9 and coef<1e-9 and norm<1e-9))
        results.append(item);print(json.dumps(item),flush=True)
    report['checks']=results;report['passed']=all(r['passed'] for r in results)
elif a.mode=='derivatives':
    from methanol_release.forward import DATA,load_problem,reduced_transform
    endpoint=np.load(DATA/'endpoint.npz',allow_pickle=False);x=endpoint['x']
    problem=load_problem();values=problem.evaluate(x,True)
    J=values['jacobian'].copy();T=reduced_transform()
    # All independent molecular directions plus mixed response/full directions.
    rng=np.random.default_rng(20261108);directions=[]
    norms=np.maximum(np.linalg.norm(J@T,axis=0),1.)
    scaled=T/norms
    for i in range(13):directions.append((f'molecular-{i}',scaled[:,i]))
    for i in range(6):
        z=rng.normal(size=127);z[:13]=0
        directions.append((f'response-{i}',scaled@z/np.linalg.norm(z)))
    for i in range(3):
        z=rng.normal(size=127)
        directions.append((f'mixed-{i}',scaled@z/np.linalg.norm(z)))
    checks=[]
    shape_indices=[problem.global_slice.start+s['slice'].start+2
                   for s in problem.global_problem.specs]
    for name,v in directions:
        expected=J@v
        radius=min([min(x[i],1-x[i])/abs(v[i]) for i in shape_indices if v[i]] or [np.inf])
        # Keep both evaluations inside the physical Voigt-shape domain.
        # Molecular bounds constrain optimization, not the derivative domain.
        step=min(1e-3,.25*radius)
        if step<=1e-10:raise ArithmeticError('No resolvable two-sided derivative step')
        for h in [step,step/2]:
            fd=(problem.residual(x+h*v)-problem.residual(x-h*v))/(2*h)
            error=float(np.max(abs(fd-expected))/max(1.,np.max(abs(expected))))
            checks.append(dict(direction=name,step=h,max_scaled_error=error,passed=error<5e-5))
    report.update(checks=checks,comparison_max=float(np.max(abs(J-endpoint['jacobian']))),
        selective_table_comparison_max=float(np.max(abs(J[:99]-endpoint['jacobian'][:99]))),
        max_independent_difference=max(r['max_scaled_error'] for r in checks),
        threshold=5e-5,passed=all(r['passed'] for r in checks))
else:
    from methanol_release.bootstrap import QuadraticBootstrap,noise_stream
    model=QuadraticBootstrap();checks=[]
    for label,noise in [('zero',np.zeros(10142)),('nominal-draw-0',next(noise_stream())[1])]:
        status,delta=model.solve(noise)
        checks.append(dict(label=label,**status))
    report.update(linearization=model.wavefunction_source,allowance=model.allowance,
        checks=checks,passed=all(c['valid'] for c in checks))
report['seconds']=time.monotonic()-start
a.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('configuration','checks')},indent=2))
raise SystemExit(0 if report['passed'] else 1)
