#!/usr/bin/env python3
"""Recalculate and render all L26 paper figures and tables from public inputs."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--main-python', default=sys.executable,
                        help='Python with NumPy 2.3.5 and Matplotlib 3.11.2')
    parser.add_argument('--atlas-python', default=sys.executable,
                        help='Python with Matplotlib 3.10.7 (the numerical environment)')
    args = parser.parse_args()
    out = args.output.absolute()
    if any(p.is_symlink() for p in [out, *out.parents]):
        parser.error('Symlink output refused; use the real absolute path')
    if out.exists() or out.resolve().is_relative_to(ROOT):
        parser.error('Choose a new output directory outside the source tree')
    out.mkdir(parents=True)
    numerical = out/'numerical'
    os.environ.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
        OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONPATH=str(ROOT/'src'),
        METHANOL_BELGI_CACHE=str(numerical/'native-cache'), MPLBACKEND='Agg')
    sys.path.insert(0,str(ROOT/'src'))
    calls = []

    def run(label, command):
        print(label, flush=True)
        with (out/(label+'.log')).open('w') as log:
            result = subprocess.run(command,cwd=ROOT,env=os.environ,text=True,stdout=log,stderr=subprocess.STDOUT)
        calls.append(dict(step=label,command=command,exit_code=result.returncode))
        if result.returncode:
            raise RuntimeError('Reproduction failed; see '+str(out/(label+'.log')))

    run('numerical', [sys.executable,str(ROOT/'reproduce.py'),'--output',str(numerical)])
    run('uncertainties', [sys.executable,'-m','methanol_release.bootstrap','--run',
                          '--output',str(out/'uncertainties')])
    import numpy as np
    from methanol_release.forward import DATA, load_problem
    from paper_tables import build_tables
    from paper_figure_data import build_figure_data
    problem = load_problem()
    with np.load(DATA/'endpoint.npz',allow_pickle=False) as archive:
        coordinates = archive['x'].copy()
    print('tables',flush=True)
    tables = build_tables(out/'tables',problem,coordinates,out/'uncertainties/draws.npz',numerical)
    print('figure-calculations',flush=True)
    figure_source = out/'figure-source/current-figures'
    calculations = build_figure_data(figure_source,problem,coordinates,numerical)
    run('render', [sys.executable,str(figure_source/'reproduce.py'),'--output',str(out/'figures'),
                   '--main-python',args.main_python,'--atlas-python',args.atlas_python])
    figures = []
    for path in sorted((out/'figures').glob('*/figures/*.pdf')):
        relative = path.relative_to(out/'figures')
        reference = ROOT/'current-figures'/relative
        figures.append(dict(path=str(path.relative_to(out)),sha256=sha256(path),
            identical_to_publication=sha256(path)==sha256(reference)))
    report = dict(passed=tables['passed'] and calculations['passed'] and len(figures)==12,
        model='L26',scope='Recompute published-model observables, conditional draws and rendering; no molecular refit',
        source_manifest_sha256=sha256(ROOT/'MANIFEST.json'),
        python=platform.python_version(),platform=platform.platform(),
        packages={name:importlib.metadata.version(name) for name in ['numpy','scipy','sympy','matplotlib','pywigxjpf']},
        tables=tables,figure_calculations=calculations,figures=figures,calls=calls,
        comparison_policy='Exact table text and numerical tolerances are required. PDF byte matches are reported separately because renderers/fonts can differ.')
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['passed']:
        raise RuntimeError('Publication comparison failed; see '+str(out/'validation.json'))
    print('Reproduced 12 figures, 8 main tables and 27 supplementary tables in '+str(out))


if __name__=='__main__':
    main()
