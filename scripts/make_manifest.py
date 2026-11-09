#!/usr/bin/env python3
"""Seal the explicitly reviewed inventory after numerical and export checks."""
import json,hashlib
from pathlib import Path
from export_safety import selected_paths,safe_file,digest
ROOT=Path(__file__).resolve().parents[1]
rows=[dict(path=n,bytes=safe_file(ROOT,n).stat().st_size,sha256=digest(safe_file(ROOT,n))) for n in selected_paths(ROOT)]
(ROOT/'MANIFEST.json').write_text(json.dumps(dict(fit='L26',files=rows),indent=2)+'\n')
(ROOT/'SHA256SUMS').write_text(''.join(r['sha256']+'  '+r['path']+'\n' for r in rows))
print('Sealed',len(rows),'reviewed files')
