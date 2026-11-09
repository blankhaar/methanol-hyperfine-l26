from pathlib import Path
import csv, hashlib, json, platform
from collections import defaultdict
import numpy as np
import scipy
from .profile_estimators import SHAPES,fitted_center,broad_limit,component_width_for_composite,component_width_for_fitted_fwhm,absorption_fwhm
HERE=Path(__file__).resolve().parent/'data'
from .benchmarks import BENCHMARKS,L1498_WIDTHS
TARGETS={r['target']:r['id'] for r in BENCHMARKS}
FREQ_GHZ={r['target']:r['frequency_MHz']/1000 for r in BENCHMARKS}
# Illustrative per-component widths. No claim that one equals any observed scan.
WIDTHS=(30.,100.,300.)

def generate(output):
    global HERE
    HERE=Path(output)
    groups=defaultdict(list)
    for r in csv.DictReader((HERE/'components.csv').open()):groups[r['fit'],r['transition']].append(r)
    rows=[];tests=[]
    for (fit,tr),rr in groups.items():
        if tr not in TARGETS.values():continue
        target=next(k for k,v in TARGETS.items() if v==tr)
        f=np.array([float(r['offset_kHz']) for r in rr]);s=np.array([float(r['strength']) for r in rr]);s/=s.sum()
        for shape in SHAPES:
            for w in WIDTHS:
                r=fitted_center(f,s,w,shape)
                rows.append(dict(fit=fit,target=target,transition=tr,response=shape,component_absorption_FWHM_kHz=w,composite_absorption_FWHM_kHz=absorption_fwhm(f,s,w,shape),fitted_center_Hz=r['center']*1000,fitted_absorption_FWHM_kHz=r['fitted_absorption_FWHM'],broad_limit_Hz=broad_limit(f,s,w,shape)*1000,stationarity_norm=r['stationarity_norm'],scenario='controlled component width'))
        # At 300 kHz composite absorption width, compare the response operators fairly.
        for shape in SHAPES:
            w=component_width_for_composite(f,s,300.,shape)
            r=fitted_center(f,s,w,shape)
            rows.append(dict(fit=fit,target=target,transition=tr,response=shape,component_absorption_FWHM_kHz=w,composite_absorption_FWHM_kHz=300.,fitted_center_Hz=r['center']*1000,fitted_absorption_FWHM_kHz=r['fitted_absorption_FWHM'],broad_limit_Hz=broad_limit(f,s,w,shape)*1000,stationarity_norm=r['stationarity_norm'],scenario='controlled composite absorption width'))
        if target in ('12.2GHz','48.372GHz'):
            # Rest-frame equivalent of the published 12 km/s FWHM; all components
            # assigned statistical populations; component width is chosen so the
            # single-Gaussian fitted FWHM matches the quoted width scale.
            gamma=FREQ_GHZ[target]*1e6*12/299792.458
            w=component_width_for_fitted_fwhm(f,s,gamma,'Gaussian')
            r=fitted_center(f,s,w,'Gaussian')
            rows.append(dict(fit=fit,target=target,transition=tr,response='Gaussian',component_absorption_FWHM_kHz=w,composite_absorption_FWHM_kHz=absorption_fwhm(f,s,w,'Gaussian'),fitted_center_Hz=r['center']*1000,fitted_absorption_FWHM_kHz=r['fitted_absorption_FWHM'],broad_limit_Hz=broad_limit(f,s,w,'Gaussian')*1000,stationarity_norm=r['stationarity_norm'],scenario='PKS1830 fitted FWHM 12 km per s'))
        if target in L1498_WIDTHS:
            gamma=L1498_WIDTHS[target]
            w=component_width_for_fitted_fwhm(f,s,gamma,'Gaussian')
            r=fitted_center(f,s,w,'Gaussian')
            rows.append(dict(fit=fit,target=target,transition=tr,response='Gaussian',component_absorption_FWHM_kHz=w,composite_absorption_FWHM_kHz=absorption_fwhm(f,s,w,'Gaussian'),fitted_center_Hz=r['center']*1000,fitted_absorption_FWHM_kHz=r['fitted_absorption_FWHM'],broad_limit_Hz=broad_limit(f,s,w,'Gaussian')*1000,stationarity_norm=r['stationarity_norm'],scenario='L1498-1 published fitted width scale'))
    with (HERE/'profile_results.csv').open('w') as o:
        wr=csv.DictWriter(o,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    (HERE/'profile_results.json').write_text(json.dumps(dict(rows=rows,units='kHz internally; center Hz',python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,source_catalog_sha256=hashlib.sha256((HERE/'components.csv').read_bytes()).hexdigest()),indent=2)+'\n')
    print(f'Computed {len(rows)} L26 response scenarios.')
