"""Harmless synthetic fixtures exercising the export trust boundary."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from export_safety import validate_names,private_reason,safe_file,scan,verify,digest,selected_paths

class ExportSecurity(unittest.TestCase):
 def test_unsafe_paths(self):
   for names in [['/tmp/escape'],['../escape'],['a/../b'],['a//b'],['./a'],['a\\b'],['C:/escape'],['a','a']]:
    with self.subTest(names=names),self.assertRaises(ValueError):validate_names(names)
 def test_private_categories(self):
   for name in ['publication/paper.pdf','private-audit/a.json','private-source-data/a.csv','validation/a.json','supplement/tex/a.pdf','a/.env','a/.env.local','draft-manuscript.tex','release-evidence/claude-access.json']:
    self.assertTrue(private_reason(name),name)
 def test_symlink_parent(self):
   with tempfile.TemporaryDirectory() as d:
    root=Path(d).resolve();(root/'real').mkdir();(root/'real/a').write_text('harmless');(root/'link').symlink_to(root/'real',target_is_directory=True)
    with self.assertRaises(ValueError):safe_file(root,'link/a')
 def test_large_file_scanned(self):
   with tempfile.TemporaryDirectory() as d:
    p=Path(d).resolve()/'fixture';p.write_bytes(b'x'*3_000_000+b'\n'+b'-----BEGIN '+b'PRIVATE KEY-----')
    self.assertIn('private_key',scan(p))
 def test_unknown_files_not_selected_and_export_rejects(self):
   with tempfile.TemporaryDirectory() as d:
    root=Path(d).resolve();(root/'export-inventory.json').write_text(json.dumps({'files':['export-inventory.json','safe.txt']}));(root/'safe.txt').write_text('harmless')
    rows=[{'path':n,'bytes':(root/n).stat().st_size,'sha256':digest(root/n)} for n in selected_paths(root)]
    (root/'MANIFEST.json').write_text(json.dumps({'files':rows}));verify(root)
    (root/'draft-manuscript.md').write_text('harmless fixture')
    self.assertNotIn('draft-manuscript.md',selected_paths(root))
    with self.assertRaises(ValueError):verify(root)
    verify(root,exact=False)
 def test_duplicate_manifest_rejected(self):
   with tempfile.TemporaryDirectory() as d:
    root=Path(d).resolve();(root/'MANIFEST.json').write_text(json.dumps({'files':[{'path':'x'},{'path':'x'}]}))
    with self.assertRaises(ValueError):verify(root)
 def test_export_rejects_before_writing_and_removes_failed_stage(self):
   from export_code import export
   from unittest.mock import patch
   with tempfile.TemporaryDirectory() as d:
    base=Path(d).resolve();root=base/'source';root.mkdir();out=base/'export'
    (root/'export-inventory.json').write_text(json.dumps({'files':['export-inventory.json','safe.txt']}));(root/'safe.txt').write_text('harmless')
    rows=[{'path':n,'bytes':(root/n).stat().st_size,'sha256':digest(root/n)} for n in selected_paths(root)]
    (root/'MANIFEST.json').write_text(json.dumps({'files':rows}));(root/'SHA256SUMS').write_text('fixture')
    with patch('export_code.audit',side_effect=[{'passed':True},{'passed':False}]):
     with self.assertRaises(ValueError):export(root,out)
    self.assertFalse(out.exists());self.assertEqual(list(base.iterdir()),[root])
    rows[0]['path']='../escape';(root/'MANIFEST.json').write_text(json.dumps({'files':rows}))
    with self.assertRaises(ValueError):export(root,out)
    self.assertFalse(out.exists());self.assertEqual(list(base.iterdir()),[root])
