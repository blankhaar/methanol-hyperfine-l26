"""Shared fail-closed export boundary. No distribution or rights authorization."""
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat

CONTROL = {'MANIFEST.json', 'SHA256SUMS'}
PRIVATE_PREFIXES = ('publication/', 'private-source-data/', 'private-audit/', 'local-only-upstream/', 'validation/', 'supplement/tex/', 'supplement/reference/', 'release-evidence/publication/')
PRIVATE_NAMES = {'.env', '.credentials.json', 'id_rsa', 'id_ed25519', 'claude-access.json', 'private-rights-review.json', 'private-export-audit.json'}
SKIP = {'.git', '__pycache__', '.pytest_cache', '.generated', '.venv', 'build', 'dist'}
ROOT_PATTERN = re.compile(rb'/(?:Users|home)/')
WORKSPACE_PATTERN = re.compile(rb'/' + rb'workspace/')
PATTERNS = {
 'private_key': re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
 'github_token': re.compile(rb'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{60,})\b'),
 'aws_access_key': re.compile(rb'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
 'anthropic_key': re.compile(rb'\bsk-ant-[A-Za-z0-9_-]{50,}\b'),
}

def validate_names(names):
 seen = set()
 for name in names:
  if not isinstance(name, str) or not name or '\\' in name or '\x00' in name:
   raise ValueError('Invalid path')
  p = PurePosixPath(name)
  if p.is_absolute() or name != p.as_posix() or any(x in ('', '.', '..') for x in name.split('/')) or ':' in name:
   raise ValueError('Unsafe path: ' + name)
  if name in seen: raise ValueError('Duplicate path: ' + name)
  seen.add(name)
 return seen

def private_reason(name):
 p = PurePosixPath(name)
 if name.startswith(PRIVATE_PREFIXES) or any(x in PRIVATE_NAMES or x.startswith('.env.') for x in p.parts): return 'Private category'
 if p.suffix.lower()=='.tex' and not name.startswith('tables/'): return 'Manuscript text outside export'
 if p.suffix.lower() in { '.zip', '.whl', '.pyc', '.docx', '.bundle', '.pkl', '.pickle', '.so', '.dylib', '.dll', '.exe'}: return 'Outside source-export format'
 if any(x in SKIP or x.endswith('.egg-info') for x in p.parts): return 'Build/cache/history category'
 return None

def safe_file(root, name):
 validate_names([name]); root = Path(root).absolute(); p = root
 for ancestor in (root, *root.parents):
  if ancestor.is_symlink(): raise ValueError('Symlink root ancestry refused')
 for part in PurePosixPath(name).parts:
  p = p / part
  if p.is_symlink(): raise ValueError('Symlink refused: ' + name)
 if not p.is_file() or not stat.S_ISREG(p.stat().st_mode): raise ValueError('Nonregular/missing file: ' + name)
 if not p.resolve().is_relative_to(root.resolve()): raise ValueError('Escaping file: ' + name)
 return p

def selected_paths(root):
 names = json.loads(safe_file(root, 'export-inventory.json').read_text())['files']
 validate_names(names)
 for name in names:
  reason = private_reason(name)
  if reason: raise ValueError(reason + ': ' + name)
 return sorted(names)

def digest(path):
 h = hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
 return h.hexdigest()

def scan(path):
 """Scan every byte in bounded memory, including matches crossing blocks."""
 found = set(); tail = b''
 with path.open('rb') as f:
  for block in iter(lambda: f.read(1024 * 1024), b''):
   data = tail + block
   for label, pattern in PATTERNS.items():
    if pattern.search(data): found.add(label)
   if ROOT_PATTERN.search(data): found.add('personal_absolute_root')
   if WORKSPACE_PATTERN.search(data): found.add('historical_workspace_root')
   tail = data[-4096:]
 return sorted(found)

def tree_files(root, private=False):
 names = []
 for directory, dirs, files in os.walk(root, followlinks=False):
  for n in dirs + files:
   p = Path(directory) / n
   if p.is_symlink(): raise ValueError('Symlink refused: ' + p.relative_to(root).as_posix())
  if private: dirs[:] = [n for n in dirs if n not in SKIP and not n.endswith('.egg-info')]
  for n in files:
   if private and n == '.DS_Store': continue
   names.append((Path(directory) / n).relative_to(root).as_posix())
 return sorted(names)

def verify(root, manifest_name='MANIFEST.json', exact=True):
 manifest = json.loads(safe_file(root, manifest_name).read_text()); rows = manifest['files']
 names = [r['path'] for r in rows]; validate_names(names)
 if manifest_name == 'MANIFEST.json':
  if set(names) != set(selected_paths(root)): raise ValueError('Manifest differs from reviewed inventory')
  for name in names:
   if private_reason(name): raise ValueError('Private manifest member: ' + name)
 for r in rows:
  p = safe_file(root, r['path'])
  if p.stat().st_size != r['bytes'] or digest(p) != r['sha256']: raise ValueError('Size/hash mismatch: ' + r['path'])
 if exact:
  extras = set(tree_files(root)) - set(names) - CONTROL
  if extras: raise ValueError('Unlisted files: ' + ', '.join(sorted(extras)))
 return rows
