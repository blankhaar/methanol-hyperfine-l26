#!/usr/bin/env python3
"""Create an atomic fresh local candidate. No upload or push is performed."""
import argparse,json,os,shutil,subprocess,tempfile
from pathlib import Path
import sys
sys.dont_write_bytecode = True
from export_safety import verify, safe_file
from audit_private_export import audit
ROOT=Path(__file__).resolve().parents[1]

def export(root,out):
 root=Path(root).resolve();out=Path(out).absolute()
 # Reject every symlink component before resolving or writing anything.
 for p in [out,*out.parents]:
  if p.is_symlink():raise ValueError('Output symlink refused')
 if out.exists() or out.resolve().is_relative_to(root):raise ValueError('Output must be a fresh path outside source')
 if not out.parent.is_dir():raise ValueError('Output parent must already exist')
 rows=verify(root,exact=False)
 report=audit(root,export=False)
 if not report['passed']:raise ValueError('Source privacy preflight failed: '+json.dumps(report))
 staging=Path(tempfile.mkdtemp(prefix='.'+out.name+'-',dir=out.parent))
 reserved=False
 try:
  for name in [r['path'] for r in rows]+['MANIFEST.json','SHA256SUMS']:
   target=staging/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(safe_file(root,name),target)
  verify(staging)
  report=audit(staging,export=True)
  if not report['passed']:raise ValueError('Staged privacy audit failed: '+json.dumps(report))
  # Execute the actual detached payload before publishing the destination.
  with tempfile.TemporaryDirectory(prefix='methanol-export-tests-') as test_work:
   env=os.environ.copy();env.pop('PYTHONPATH',None);env.pop('METHANOL_TEST_PRIVATE_DOCUMENTS',None)
   env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg',MPLCONFIGDIR=test_work)
   tested=subprocess.run([sys.executable,str(staging/'scripts/run_tests.py')],cwd=test_work,env=env,text=True,capture_output=True,timeout=300)
   if tested.returncode:
    raise ValueError('Detached export tests failed:\n'+tested.stdout+tested.stderr)
   report['detached_tests']={'passed':True,'output':tested.stdout+tested.stderr}
  verify(staging)  # Tests must not change the selected payload or add files.
  # mkdir reserves the final name exclusively; rename replaces only our empty reservation.
  out.mkdir();reserved=True;os.replace(staging,out);reserved=False
 except BaseException:
  shutil.rmtree(staging,ignore_errors=True)
  if reserved:out.rmdir()
  raise
 return report

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 try:export(ROOT,a.output)
 except (ValueError,OSError,subprocess.TimeoutExpired) as e:p.exit(1,str(e)+'\n')
 print('L26 source export verified; no upload or push performed.')
if __name__=='__main__':main()
