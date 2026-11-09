"""Native compilation, eigenfunction provenance and downstream L26 checks."""
from pathlib import Path
import json
import os
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from methanol_hyperfine.belgi import PROJECT, State, export_wavefunctions
from methanol_hyperfine.native_inputs import HamiltonianError


class NativeBelgi(unittest.TestCase):
    def test_fresh_native_run_without_reading_packaged_waves(self):
        states = [State('A',0,5,1,'+'), State('A',0,6,0,'+'),
                  State('E',0,5,1), State('E',1,6,-2)]
        reference = export_wavefunctions(states, backend='reference')
        real_load = np.load
        def guard(path, *args, **kwargs):
            if isinstance(path, (str,Path)) and Path(path).resolve().is_relative_to(PROJECT/'data/wavefunctions'):
                raise AssertionError('Native route attempted to read packaged eigenfunctions')
            return real_load(path,*args,**kwargs)
        with tempfile.TemporaryDirectory() as work, patch('numpy.load',side_effect=guard):
            work = Path(work).resolve()
            native = export_wavefunctions(states,force=True,cache_dir=work)
            cached = export_wavefunctions(states,cache_dir=work)
            manifests = list(Path(work).glob('*/waves-k10-*.json'))
            self.assertEqual(len(manifests),1)
            meta = json.loads(manifests[0].read_text())
            self.assertFalse(meta['archived_wavefunctions_used'])
            self.assertEqual(meta['source_files'],58)
            self.assertIn('setup.f',meta['source_patches'])
            for state in states:
                self.assertLess(abs(native[state].energy-reference[state].energy),1e-9)
                a,b=native[state].coefficients,reference[state].coefficients
                self.assertLess(min(np.max(abs(a-b)),np.max(abs(a+b))),1e-9)
                np.testing.assert_array_equal(native[state].coefficients,cached[state].coefficients)

    def test_missing_compiler_does_not_fall_back(self):
        with patch.dict(os.environ,{'METHANOL_FORTRAN_COMPILER':'methanol-no-such-compiler'}):
            with self.assertRaisesRegex(HamiltonianError,'Archived waves are not used'):
                export_wavefunctions([State('A',0,0,0,'+')])

    def test_native_source_change_is_detected_before_execution(self):
        from methanol_hyperfine import native_belgi
        with tempfile.TemporaryDirectory() as work:
            path=Path(work);(path/'source').mkdir();(path/'source/setup.f').write_text('invalid fixture')
            with patch.object(native_belgi,'DATA',path):
                with self.assertRaisesRegex(HamiltonianError,'source differs'):
                    native_belgi.configuration()

    def test_native_fit_rebuild_matches_endpoint(self):
        from methanol_release.fit import audit
        report=audit(True)
        self.assertEqual(report['wavefunction_provider'],'native')
        self.assertLess(report['max_residual_error'],1e-10)
        self.assertAlmostEqual(report['Q_total'],2521.2262069897515,places=8)

    def test_native_likelihood_replaces_saved_operators(self):
        from methanol_release.forward import load_problem
        p=load_problem()
        self.assertEqual(p.d.wavefunction_source,'native')
        self.assertEqual(len(p.d.prepared),204)

    def test_catalog_and_line_shifts_do_not_open_saved_waves(self):
        from methanol_zeeman.catalog import _inputs,compute_transition
        from methanol_line_shift.estimate import _native_components,load_components
        _inputs.cache_clear();_native_components.cache_clear()
        real_load=np.load
        prohibited={'available_wavefunctions.npz','maser_frozen_states.npz',
                    'reference_hamiltonians.npz','complete_frozen_projections.npz'}
        def guard(path,*args,**kwargs):
            if isinstance(path,(str,Path)):
                q=Path(path).resolve()
                if q.name in prohibited or q.is_relative_to(PROJECT/'data/wavefunctions'):
                    raise AssertionError('Native calculation attempted a reference wave/H read')
            return real_load(path,*args,**kwargs)
        with patch('numpy.load',side_effect=guard):
            self.assertTrue(compute_transition('6.7'))
            offsets,strengths=load_components('12.2GHz')
            self.assertTrue(np.isfinite(offsets).all())
            self.assertGreater(strengths.sum(),0.)


if __name__=='__main__':
    unittest.main()
