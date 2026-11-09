"""Reproduce L26, figures, tensors, magnetic coefficients and line shifts."""
import argparse,json,subprocess,sys
from pathlib import Path

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['belgi','fit','bootstrap','zeeman','shift','hamiltonian','figures','verify','operators','status'])
    a,rest=p.parse_known_args(argv)
    if a.command == 'belgi':
        from methanol_hyperfine.belgi import main as run
        return run(rest)
    if a.command in ['fit','bootstrap']:
        from importlib import import_module
        return import_module('methanol_release.'+a.command).main(rest)
    if a.command=='zeeman':
        from methanol_zeeman.__main__ import main as run
        return run(rest)
    if a.command=='shift':
        from methanol_line_shift.__main__ import main as run
        return run(rest)
    if a.command=='operators':
        if rest:p.error('operators takes no arguments')
        from .operators import verify_operators
        print(json.dumps(verify_operators(),indent=2));return 0
    if a.command=='hamiltonian':
        q=argparse.ArgumentParser();q.add_argument('--state',default='A_0_5_1_+');q.add_argument('--output',type=Path,required=True);x=q.parse_args(rest)
        from methanol_zeeman.catalog import _inputs
        from methanol_zeeman.frozen_model import hyperfine,parameters
        from methanol_zeeman.output_paths import validate_output_file
        import numpy as np
        output=validate_output_file(x.output)
        if output.suffix!='.npz':q.error('Output must end in .npz')
        waves=_inputs()['available_wavefunctions']
        if x.state not in waves:q.error('State unavailable')
        h=hyperfine(x.state,waves[x.state]);output.parent.mkdir(parents=True,exist_ok=True)
        with output.open('xb') as f:np.savez_compressed(f,hamiltonian_kHz=h,state=x.state,fit_sha256=parameters['source_sha256'])
        return 0
    root=Path(__file__).resolve().parents[2]
    if a.command=='status':
        from importlib.resources import files
        print(files('methanol_release').joinpath('data/package-status.txt').read_text(),end='');return 0
    target='current-figures/reproduce.py' if a.command=='figures' else 'scripts/verify_manifest.py'
    return subprocess.call([sys.executable,str(root/target),*rest])
if __name__=='__main__':raise SystemExit(main())
