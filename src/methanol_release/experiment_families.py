"""Recompute Table VIII with an explicit positive-response peak-selection rule."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar
from .forward import DATA
from methanol_hyperfine.lambdip import lorentz_2f
from methanol_hyperfine.belgi import State, export_wavefunctions
from methanol_zeeman.frozen_model import hyperfine
from methanol_zeeman.spin import transition_components
from methanol_zeeman.output_paths import validate_output_directory

def peaks(f,s,step=.02,height=.10,prominence=.05):
    """Positive central 2f phase; thresholds relative to the profile maximum."""
    x=np.arange(float(f.min())-50,float(f.max())+50,step)
    y=-lorentz_2f(x[:,None]-f,5.,3.)@s
    maximum=float(y.max())
    ids,_=find_peaks(y,height=height*maximum,prominence=prominence*maximum)
    def response(v):return -float(lorentz_2f(v-f,5.,3.)@s)
    positions=np.array([minimize_scalar(lambda v:-response(v),bounds=(x[i-1],x[i+1]),method='bounded',options={'xatol':1e-10}).x for i in ids])
    return positions

def calculate(output):
    out=validate_output_directory(output,('transitions.json','families.csv','validation.json'));out.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((DATA/'experiment_families.csv').open()))
    states = {State.from_key(r[q]) for r in rows for q in ['upper','lower']}
    waves = export_wavefunctions(states)
    matrices = {s.key: hyperfine(s.key, w.coefficients) for s,w in waves.items()}
    result=[];grid_error=0.
    for row in rows:
        u,l=row['upper'],row['lower'];c=transition_components(int(u.split('_')[2]),int(l.split('_')[2]),u[0],matrices[u],matrices[l])
        f,s=c[:,:2].T;s=s/s.sum();p=peaks(f,s);fine=peaks(f,s,step=.01)
        if len(p)!=len(fine):raise ArithmeticError('Peak count changes with sampling grid')
        error=float(np.max(abs(p-fine))) if len(p) else 0.;grid_error=max(grid_error,error)
        if error>1e-4:raise ArithmeticError('Peak position not converged')
        sensitivity={}
        for tag,h,pr in [('loose',.05,.025),('strict',.20,.10)]:
            q=peaks(f,s,height=h,prominence=pr);sensitivity[tag]=dict(count=len(q),gaps_kHz=np.diff(q).tolist())
        result.append(dict(**row,count=len(p),peaks_kHz=p.tolist(),gaps_kHz=np.diff(p).tolist(),threshold_sensitivity=sensitivity))
    families=[]
    for tag in dict.fromkeys(r['tag'] for r in result):
        rr=[r for r in result if r['tag']==tag];g=[x for r in rr for x in r['gaps_kHz']]
        if not g:raise ArithmeticError('Family has no resolved interval')
        families.append(dict(row=int(rr[0]['row']),family=tag,transitions=len(rr),N_min=min(r['count'] for r in rr),N_max=max(r['count'] for r in rr),gap_min_kHz=min(g),gap_max_kHz=max(g),single_feature_lower_J=','.join(r['new_lower_state_J'] for r in rr if r['count']==1)))
    (out/'transitions.json').write_text(json.dumps(result,indent=2)+'\n')
    with (out/'families.csv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(families[0]),lineterminator='\n');w.writeheader();w.writerows(families)
    report=dict(passed=True,fit='L26',wavefunction_provider='native BELGI',families=len(families),transitions=len(rows),maximum_grid_position_error_kHz=grid_error,
        rule=dict(Lorentz_HWHM_kHz=5.,modulation_amplitude_kHz=3.,phase='positive isolated-line center',height_fraction=.10,prominence_fraction=.05,grid_kHz=.02,verification_grid_kHz=.01,statistical_strengths=True,crossovers=False),
        scope='Declared numerical feature-selection rule; not calibrated experimental resolution',threshold_sensitivity='Separate 5%/2.5% and 20%/10% height/prominence results retained for every transition')
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    print(json.dumps(calculate(p.parse_args().output),indent=2))
