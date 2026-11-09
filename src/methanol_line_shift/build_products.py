"""Recompute response scenarios and bounded rational line-pull approximations."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .rebuild import rebuild_response_catalog
from .profile_estimators import SHAPES,broad_limit,fitted_center
from .run_profile_audit import generate as scenarios
from .mass_ratio_scenarios import generate as mass_ratio
from methanol_zeeman.output_paths import validate_output_directory
DEGREES=[(0,1),(1,1),(1,2),(2,2),(2,3),(3,3),(3,4),(4,4),(4,5),(5,5),(5,6),(6,6)]
def rounded(x):
    return float(f"{float(x):.12g}")

def evaluate(w, model):
    x = (30. / w)**2
    return (model["A_Hz_kHz2"] / w**2 *
            np.polynomial.polynomial.polyval(x, model["numerator"]) /
            np.polynomial.polynomial.polyval(x, model["denominator"]))

def poly_min(coefficients):
    """Global polynomial minimum on 0 <= x <= 1, including stationary points."""
    roots = np.polynomial.polynomial.polyroots(
        np.polynomial.polynomial.polyder(coefficients))
    candidates = [0., 1.] + [float(z.real) for z in roots
                             if abs(z.imag) < 1e-10 and 0 < z.real < 1]
    return float(np.min(np.polynomial.polynomial.polyval(candidates, coefficients)))

def fit_rational(w, y, a, m, n):
    x = (30. / w)**2
    r = y * w**2 / a
    X = np.array([x**j for j in range(1, m+1)]).T if m else np.empty((len(x), 0))
    Z = np.array([x**j for j in range(1, n+1)]).T
    initial = np.linalg.lstsq(np.column_stack((X, -r[:, None]*Z)), r-1, rcond=None)[0]

    def residual(z):
        return ((1+X@z[:m]) / (1+Z@z[m:]) - r) / r

    result = least_squares(residual, initial, xtol=1e-14, ftol=1e-14,
                           gtol=1e-14, max_nfev=1000)
    model = dict(A_Hz_kHz2=rounded(a), numerator=[1.] + list(map(rounded, result.x[:m])),
                 denominator=[1.] + list(map(rounded, result.x[m:])),
                 numerator_degree=m, denominator_degree=n,
                 coefficient_significant_digits=12)
    model["minimum_denominator_x_0_to_1"] = poly_min(model["denominator"])
    model["minimum_numerator_x_0_to_1"] = poly_min(model["numerator"])
    return model

def errors(w, exact, model):
    residual = evaluate(w, model) - exact
    rel = np.abs(residual / exact)
    i, j = np.argmax(np.abs(residual)), np.argmax(rel)
    return dict(max_abs_error_Hz=float(abs(residual[i])),
                width_at_max_abs_error_kHz=float(w[i]),
                max_relative_error=float(rel[j]),
                width_at_max_relative_error_kHz=float(w[j]),
                rms_error_Hz=float(np.sqrt(np.mean(residual**2))))

def generate(output):
    out=validate_output_directory(output,('components.csv','profile_results.csv','line_pull_surrogate_coefficients.json'))
    rebuild_response_catalog(out)
    scenarios(out);mass_ratio(out)
    components=list(csv.DictReader((out/'components.csv').open()));models={}
    for target,transition in [('12.2GHz','C2015-p06-1'),('96.741GHz','benchmark-96.741GHz')]:
        rows=[r for r in components if r['transition']==transition]
        f=np.array([float(r['offset_kHz']) for r in rows]);s=np.array([float(r['strength']) for r in rows]);s/=s.sum()
        widths=np.geomspace(30,2000,183);validation=np.geomspace(30.017,1999.13,151)
        for shape in SHAPES:
            y=np.array([1000*(fitted_center(f,s,w,shape)['center']-s@f) for w in widths])
            exact=np.array([1000*(fitted_center(f,s,w,shape)['center']-s@f) for w in validation])
            leading=1000*broad_limit(f,s,1.,shape)
            for m,n in DEGREES:
                model=fit_rational(widths,y,leading,m,n);test=errors(validation,exact,model)
                if model['minimum_denominator_x_0_to_1']>0 and test['max_abs_error_Hz']<1e-3 and test['max_relative_error']<1e-6:break
            else:raise ArithmeticError('No validated surrogate for '+target+'/'+shape)
            model['validation']=test;models[target+'/'+shape]=model
    (out/'line_pull_surrogate_coefficients.json').write_text(json.dumps(dict(fit='L26',formula='c_Hz=(A/W**2)*P((30/W)**2)/Q((30/W)**2)',component_absorption_FWHM_domain_kHz=[30,2000],models=models,
        uncertainty_note='Approximation errors only; conditional on L26 and response assumptions'),indent=2)+'\n')
    return dict(passed=True,wavefunction_provider='native BELGI',models=len(models),scenarios=232)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    print(json.dumps(generate(p.parse_args().output),indent=2))
