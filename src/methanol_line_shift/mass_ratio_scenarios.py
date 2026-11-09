"""Conditional common-velocity projection, NOT a revised observed mu constraint.

mu=m_p/m_e, K=d ln nu/d ln mu, epsilon=(h_astro-h_ref)/nu.
Fit epsilon_i = b + K_i * delta_mu. b absorbs a common fractional Doppler shift.
The source-specific h_ref values are UNKNOWN; zero is set only to isolate the
simulated astronomical-side contribution. Quoted frequency-center errors are
used as independent observation errors; no reference uncertainty/covariance,
source kinematics, line-width uncertainty, or actual radiative transfer is added.
"""
from pathlib import Path
import csv,json
import numpy as np
from .benchmarks import BENCHMARKS,L1498_OBS_SIGMA,L1498_XU_MHZ
HERE=Path(__file__).resolve().parent/'data'
B={x['target']:x for x in BENCHMARKS}

def projected_bias(targets,h,frequency_sigma_Hz,frequency_MHz=None):
    nu=np.array([(frequency_MHz[t] if frequency_MHz else B[t]['frequency_MHz'])*1e6 for t in targets]);K=np.array([B[t]['K_mu'] for t in targets]);eps=np.array([h[t] for t in targets])/nu;sig=np.array([frequency_sigma_Hz[t] for t in targets])/nu
    # Rescale all uncertainties by a common value for numerical conditioning.
    scale=float(np.min(sig));weights=scale/sig
    X=np.c_[np.ones(len(K)),K];A=X*weights[:,None]
    beta=np.linalg.lstsq(A,eps*weights,rcond=None)[0]
    cov=np.linalg.inv(A.T@A)*scale**2
    shifted=np.linalg.lstsq(A,(eps+1e-7)*weights,rcond=None)[0]
    assert abs(shifted[1]-beta[1])<1e-20
    if len(targets)==2:
        exact=(eps[0]-eps[1])/(K[0]-K[1]);u=np.hypot(*sig)/abs(K[0]-K[1])
        assert abs(beta[1]-exact)<1e-20 and abs(np.sqrt(cov[1,1])-u)<1e-20
    return dict(targets=targets,apparent_delta_mu_over_mu=float(beta[1]),common_fractional_offset=float(beta[0]),observed_center_only_sigma_mu=float(np.sqrt(cov[1,1])),absolute_bias_over_observed_center_only_sigma=float(abs(beta[1])/np.sqrt(cov[1,1])))

def generate(output):
    global HERE
    HERE=Path(output)
    rows=list(csv.DictReader((HERE/'profile_results.csv').open()));out=[]
    spec=[('L1498-1 published fitted width scale', [['96.741GHz','157.271GHz'],['96.739GHz','157.271GHz'],['96.739GHz','96.741GHz','157.271GHz'],list(L1498_OBS_SIGMA)],L1498_OBS_SIGMA),
          ('PKS1830 fitted FWHM 12 km per s',[['12.2GHz','48.372GHz']],{t:B[t]['frequency_MHz']*1e6*v/299792.458 for t,v in [('12.2GHz',.7),('48.372GHz',.1)]})]
    for fit in sorted(set(r['fit'] for r in rows)):
        for scenario,sets,errors in spec:
            h={r['target']:float(r['fitted_center_Hz']) for r in rows if r['fit']==fit and r['scenario']==scenario}
            for targets in sets:out.append(dict(fit=fit,scenario=scenario,**projected_bias(targets,h,errors,L1498_XU_MHZ if scenario.startswith('L1498') else None)))
    report=dict(status='passed',assumptions=__doc__,published_Dapra_handling='Section6.1 already simulated broadened hyperfine spectra and Gaussian-fit pulls. They varied rest frequencies by +/-1kHz and included a maximum 1e-9 hyperfine systematic in Delta mu/mu. The present pair/three-line scenarios must not be added independently to that already-budgeted systematic or reported as a correction to their complete six-line analysis.',source='https://academic.oup.com/mnras/article/472/4/4434/4157801',rows=out,tests={'pairs_match_analytic_difference':True,'common_velocity_shift_invariance':True})
    (HERE/'mass_ratio_scenarios.json').write_text(json.dumps(report,indent=2)+'\n')
    for r in out:
        if r['fit']=='L26':print(json.dumps(r))
