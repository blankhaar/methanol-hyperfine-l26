#!/usr/bin/env python3
"""Render the five paper figures and seven supplementary plates for L26."""
import argparse,json,os,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--main-python',default=sys.executable);p.add_argument('--atlas-python',default=sys.executable)
    a=p.parse_args();out=a.output.absolute()
    if any(q.is_symlink() for q in [out,*out.parents]):raise ValueError('Symlink output refused')
    if out.exists() or out.resolve().is_relative_to(ROOT.parent):raise ValueError('Choose a new output outside source')
    out.mkdir(parents=True)
    for folder in ['scripts','evidence','inputs','supplement/scripts','supplement/data']:shutil.copytree(ROOT/folder,out/folder)
    for group in ['manuscript','supplement']:(out/group/'figures').mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MPLCONFIGDIR=str(out/'.mpl'),SOURCE_DATE_EPOCH='1793923200')
    calls=[]
    for python,script in [(a.main_python,'scripts/draw_methanol_geometry.py'),(a.main_python,'scripts/build_fit_figures.py'),(a.main_python,'scripts/build_fit_line_pulls.py'),(a.atlas_python,'supplement/scripts/build_fit_profiles.py')]:
        proc=subprocess.run([python,str(out/script)],cwd=out,env=env,capture_output=True,text=True)
        calls.append(dict(script=script,exit_code=proc.returncode))
        (out/(Path(script).stem+'.log')).write_text(proc.stdout+proc.stderr)
        if proc.returncode:raise RuntimeError(proc.stderr[-2000:])
    figures=list(out.glob('*/figures/*.pdf'))
    if len(figures)!=12:raise RuntimeError('Expected twelve figures')
    (out/'validation.json').write_text(json.dumps(dict(passed=True,fit='L26',figures=12,calls=calls),indent=2)+'\n')
    print('Rendered twelve L26 figures in '+str(out))
if __name__=='__main__':main()
