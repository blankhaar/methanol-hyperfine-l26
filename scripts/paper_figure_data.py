"""Recompute publication plotting data from the public L26 model and observations."""
from pathlib import Path
import hashlib
import json
import shutil

import numpy as np

from methanol_hyperfine.belgi import PROJECT, State, export_wavefunctions
from methanol_hyperfine.lambdip import components, prediction
from methanol_release.operators import prepare_masked_state
from methanol_line_shift.profile_estimators import SHAPES, component_width_for_fitted_fwhm
from methanol_line_shift.estimate import estimate_shift
from paper_tables import ROOT, read_csv, write_csv, response_arrays


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_figure_data(output, problem, coordinates, numerical):
    """Saved curves are comparison references only, never calculation inputs."""
    source = ROOT / 'current-figures'
    output.mkdir(parents=True)
    for folder in ['scripts', 'evidence/figure_data', 'evidence/line_pull_data',
                   'inputs', 'supplement/scripts', 'supplement/data']:
        (output/folder).mkdir(parents=True)
    # Explicitly copy only styles, labels, literature inputs and fitted parameters.
    copies = ['reproduce.py', 'inputs/fit_endpoint.json',
              'supplement/data/profile_metadata.json', 'supplement/data/endpoint_pin.json',
              'evidence/figure_data/figure_metadata.json',
              'evidence/figure_data/belov-published-table-iv.json']
    copies += [str(p.relative_to(source)) for folder in ['scripts', 'supplement/scripts']
               for p in (source/folder).glob('*.py')]
    for name in copies:
        shutil.copy2(source/name, output/name)
    data = output/'evidence/figure_data'
    atlas = output/'supplement/data'
    calculated = problem.evaluate(coordinates)
    arrays, profiles = {}, {'fit_full_weighted_residual': calculated['residual']}
    for family, records, fits in [('low', problem.d.profiles, calculated['glob']),
                                  ('direct', problem.traces, calculated['direct'])]:
        for i, (g, fit) in enumerate(zip(records, fits, strict=True)):
            key = f'{family}_{i:02}'
            mask = g.get('mask', np.ones(len(g['x']), dtype=bool))
            for suffix, value in [('x', g['x']), ('y', g['y']), ('mask', mask), ('pred', fit[1])]:
                arrays[g['id']+'__'+suffix] = value
                profiles['fit_'+key+'_prediction' if suffix=='pred' else key+'_'+suffix] = value

    hamiltonians = problem.d.hamiltonians(problem.unpack(coordinates))
    belov, _ = problem.d.belov_features(hamiltonians)
    write_csv(data/'belov.csv', [dict(ju=g['ju'], jl=g['jl'], ku=g['ku'], kl=g['kl'],
        observed_khz=g['observed_splitting_khz'], predicted_khz=float(v),
        included=bool(g['published_weight']>0)) for g,v in zip(problem.d.belov,belov,strict=True)])

    # The independent Belov triplet uses its published fitted response parameters.
    # Amplitudes/background are solved anew; no molecular optimization is performed.
    meta = json.loads((data/'figure_metadata.json').read_text())
    triplet = meta['triplet'][0]
    raw = np.genfromtxt(PROJECT/'data/processed/belov-profile-check-v1/Belov-F2.csv',
                       delimiter=',', names=True)
    g = dict(id='Belov-F2', species='E', u=State('E',0,24,1), l=State('E',0,24,0),
             x=raw['offset_khz'], y=raw['signal'], weight=triplet['weight'],
             thermal_background=True, temperature_K=293.15)
    missing = {g[q] for q in ('u','l')} - hamiltonians.keys()
    if missing:
        from methanol_hyperfine.tensors import state_hamiltonian
        for state, wave in export_wavefunctions(missing).items():
            fixed, basis = prepare_masked_state(wave, problem.d.frame)
            hamiltonians[state] = state_hamiltonian(
                (fixed, basis*problem.d.coordinate_to_fraction[:,None,None]),
                problem.unpack(coordinates))
    sticks = components(g, hamiltonians)
    residual, curve, amplitudes = prediction(g, sticks, triplet['nuisance'], True)
    arrays.update(triplet__x=g['x'], triplet__y=g['y'], triplet__with=curve)
    np.savez_compressed(data/'figure_data.npz', **arrays)
    np.savez_compressed(atlas/'fit_profiles.npz', **profiles)
    np.savez_compressed(atlas/'fit_residual.npz', residual=calculated['residual'])

    # Full least-squares response fits, with no surrogate or saved curve lookup.
    line_data = output/'evidence/line_pull_data'
    shutil.copy2(numerical/'responses/components.csv', line_data/'components.csv')
    catalog = read_csv(line_data/'components.csv')
    rows = []
    widths = np.sort(np.r_[np.geomspace(30,2000,181),100.,300.])
    for target, transition, frequency in [('12.2GHz','C2015-p06-1',12178.596415),
                                          ('96.741GHz','benchmark-96.741GHz',96741.375)]:
        f,s = response_arrays(catalog,transition)
        for shape in SHAPES:
            for width in widths:
                r = estimate_shift(f,s,width,shape)
                rows.append(dict(target=target,transition=transition,response=shape,
                    component_absorption_FWHM_kHz=float(width),
                    component_velocity_FWHM_kms=float(width*299792.458/(frequency*1000)),
                    fitted_absorption_FWHM_kHz=r['fitted_absorption_FWHM_kHz'],
                    fitted_center_minus_centroid_Hz=r['shift_Hz'],
                    stationarity_norm=r['stationarity_norm']))
    write_csv(line_data/'line_pull_curves.csv',rows)
    f,s = response_arrays(catalog,'benchmark-96.741GHz')
    width = component_width_for_fitted_fwhm(f,s,49.,'Gaussian')
    r = estimate_shift(f,s,width,'Gaussian')
    (line_data/'line_pull_figure_validation.json').write_text(json.dumps(dict(
        source_catalog_sha256=digest(line_data/'components.csv'),
        source_endpoint_sha256=meta['endpoint_sha256'],curve_rows=len(rows),
        cloud_illustration=dict(component_absorption_FWHM_kHz=width,
            fitted_blend_FWHM_kHz=r['fitted_absorption_FWHM_kHz'],
            fitted_center_minus_centroid_Hz=r['shift_Hz'],
            fitted_blend_velocity_FWHM_kms=r['fitted_absorption_FWHM_kHz']*299792.458/96741375.)),indent=2)+'\n')

    # References below are read only after every prediction has been calculated.
    with np.load(source/'evidence/figure_data/figure_data.npz',allow_pickle=False) as reference:
        differences = {k:float(np.max(abs(np.asarray(v,float)-reference[k]))) for k,v in arrays.items()}
    old_rows = read_csv(source/'evidence/line_pull_data/line_pull_curves.csv')
    def curve_key(r):
        return r['target'],r['response'],float(r['component_absorption_FWHM_kHz'])
    reference_curves = {curve_key(r):r for r in old_rows}
    line_delta = max(abs(r['fitted_center_minus_centroid_Hz']-
                        float(reference_curves[curve_key(r)]['fitted_center_minus_centroid_Hz'])) for r in rows)
    old_belov = read_csv(source/'evidence/figure_data/belov.csv')
    belov_delta = float(max(abs(v-float(r['predicted_khz'])) for v,r in zip(belov,old_belov,strict=True)))
    report = dict(passed=max(differences.values())<1e-9 and line_delta<1e-6 and belov_delta<1e-9,
        profiles=34,belov_splittings=57,line_pull_points=len(rows),
        maximum_profile_array_difference=max(differences.values()),
        maximum_belov_difference_khz=belov_delta,maximum_line_pull_difference_Hz=line_delta,
        triplet_Q=float(residual@residual),triplet_linear_amplitudes=list(map(float,amplitudes)),
        wavefunctions='fresh native BELGI', molecular_parameters='published L26 endpoint, unchanged',
        experimental_parameters='published response parameters; amplitudes/background solved again',
        reference_role='saved plotting arrays read only for comparison after calculation')
    (output/'calculation-validation.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
