#!/usr/bin/env python3
"""Check a source export for excluded documents, binaries and personal paths.

Known prohibited content identities come from the retained rights record. Path
and binary checks complement those hashes; an independent source review remains
necessary because no filename/hash list can find every unknown imported work.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import sys
sys.dont_write_bytecode = True
from export_safety import selected_paths, safe_file, private_reason, scan, tree_files, verify
ROOT = Path(__file__).resolve().parents[1]
SKIP = {'.git', '.generated', '.venv', 'build', 'dist', '__pycache__', '.pytest_cache'}


def audit(root=ROOT, export=False):
    root = Path(root).resolve()
    rights = json.loads((root / 'release-evidence/rights-audit.json').read_text())
    excluded = rights['excluded_files']
    hashes = {row['sha256']: row['path'] for row in excluded}
    approved_belgi = {r['path']: r['sha256'] for r in rights.get('bundled_belgi', {}).get('files', [])}
    forbidden_names = {Path(row['path']).name for row in excluded if 'xu2008-' in row['path']}
    matches, checked, symlinks = [], 0, []
    if export and (root / '.git').exists():
        matches.append({'path': '.git', 'reason': 'Historical Git objects must not be part of a fresh byte export'})
    try:
        names = tree_files(root) if export else selected_paths(root)
        if export:
            verify(root, exact=True)
    except (ValueError, OSError, KeyError) as e:
        return {'passed': False, 'matches': [{'reason': str(e)}], 'public_distribution_cleared': False}
    historical_paths = []; personal_paths = []
    for name in names:
        relative = Path(name)
        try:
            path = safe_file(root, name)
        except (ValueError, OSError) as e:
            matches.append({'path': name, 'reason': str(e)})
            continue
        checked += 1
        name = relative.as_posix()
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        reasons = []
        if category := private_reason(name): reasons.append(category)
        findings = scan(path)
        if 'historical_workspace_root' in findings: historical_paths.append(name)
        if 'personal_absolute_root' in findings: personal_paths.append(name)
        reasons.extend('Secret pattern: ' + x for x in findings if x not in {'personal_absolute_root','historical_workspace_root'})
        if digest in hashes:
            reasons.append('Known excluded content identity: ' + hashes[digest])
        if path.name in forbidden_names or 'vendor/belgi/' in name.lower() or '/local-cache/' in name.lower():
            reasons.append('Prohibited upstream path/name')
        if path.suffix.lower() in {'.so', '.dylib', '.dll', '.exe', '.whl', '.pkl', '.pickle', '.pyc', '.docx', '.bundle'}:
            reasons.append('Dependency binary, executable cache, publisher document or Git bundle is outside this source export')
        if name.startswith('src/methanol_hyperfine/data/belgi/') and name in approved_belgi:
            if digest != approved_belgi[name]: reasons.append('Bundled BELGI input differs from its approved identity')
        if path.suffix.lower() in {'.f', '.for', '.f90', '.f95'}:
            if approved_belgi.get(name) != digest:
                reasons.append('Fortran source is outside the exact declared BELGI inventory')
        if reasons:
            matches.append({'path': name, 'reasons': reasons})
    return {'passed': not matches and not symlinks, 'checked_files': checked,
            'known_excluded_content_hashes': len(hashes), 'matches': matches, 'symlinks': symlinks,
            'historical_absolute_path_files': historical_paths, 'personal_absolute_path_files': personal_paths, 'public_metadata_sanitization_pending': bool(personal_paths), 'export_mode': export, 'public_distribution_cleared': False,
            'scope': 'Known rights-record hashes plus explicit upstream/binary/path rules; not exhaustive rights or secret clearance'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--export', action='store_true', help='Additionally reject a .git directory')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    report = audit(args.root, args.export)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
