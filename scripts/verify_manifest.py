#!/usr/bin/env python3
"""Verify bytes, safe paths and exact export membership by default."""
import argparse,json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
from export_safety import verify

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-tree',action='store_true',help='Allow local private files outside the export inventory');p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args()
 try:
  rows=verify(a.root,exact=not a.source_tree);result={'checked':len(rows),'passed':True,'exact_export_membership':not a.source_tree}
 except (ValueError,OSError,KeyError,TypeError) as e:result={'passed':False,'error':str(e)}
 print(json.dumps(result,indent=2));return not result['passed']
if __name__=='__main__':raise SystemExit(main())
