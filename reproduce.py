#!/usr/bin/env python3
"""Validate L26 and rebuild magnetic/response outputs in a fresh directory."""
from pathlib import Path
import argparse,json,os,subprocess,sys
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
out=a.output.absolute()
if out.exists() or out.resolve().is_relative_to(ROOT):p.error('Choose a new output outside source')
if any(q.is_symlink() for q in [out,*out.parents]):p.error('Symlink output refused')
out.mkdir(parents=True)
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONPATH=str(ROOT/'src'),METHANOL_BELGI_CACHE=str(out/'native-cache'))
commands=[['scripts/verify_manifest.py','--source-tree'],['-m','methanol_hyperfine.belgi','--force','--output',str(out/'belgi')],['scripts/run_tests.py'],['-m','methanol_release.fit','--replay','--derivatives'],['-m','methanol_release.bootstrap','--replay'],['-m','methanol_zeeman','--output',str(out/'zeeman'),'--reconstruct-hyperfine'],['-m','methanol_line_shift.build_products','--output',str(out/'responses')],['-m','methanol_release.experiment_families','--output',str(out/'families')]]
results=[]
for i,cmd in enumerate(commands):
    c=subprocess.run([sys.executable,*cmd],cwd=ROOT,env=env,text=True,capture_output=True)
    (out/f'command-{i}.log').write_text(c.stdout+c.stderr);results.append(dict(command=cmd,exit_code=c.returncode))
    if c.returncode:raise RuntimeError('Failed reproduction command; see '+str(out/f'command-{i}.log'))
(out/'validation.json').write_text(json.dumps(dict(passed=True,calls=results),indent=2)+'\n')
print(out)
