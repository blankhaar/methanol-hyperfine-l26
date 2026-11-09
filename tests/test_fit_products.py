"""Scientific checks for the single published L26 endpoint."""
import csv,json,tempfile,unittest
from pathlib import Path
import numpy as np
from methanol_release.fit import audit
from methanol_release.forward import DATA,load_problem,reduced_transform
ROOT=Path(__file__).resolve().parents[1]
class FitProducts(unittest.TestCase):
    def test_exact_forward_replay_and_derivatives(self):
        report=audit(True,wavefunctions="reference")
        self.assertEqual(report['max_residual_error'],0.)
        self.assertEqual(report['max_jacobian_error'],0.)
        self.assertAlmostEqual(report['Q_total'],2521.2262069897515,places=10)
    def test_affine_basis_and_bounds(self):
        fit=json.loads((DATA/'results.json').read_text());a=np.load(DATA/'endpoint.npz',allow_pickle=False);b=np.load(DATA/'basis.npz',allow_pickle=False)
        np.testing.assert_allclose(b['M_offset']+b['M_basis']@fit['M_amplitudes'],a['x'][4:22],atol=1e-12,rtol=0)
        self.assertEqual(np.linalg.matrix_rank(reduced_transform()),127)
        self.assertTrue(np.all(a['x']>=a['lower']));self.assertTrue(np.all(a['x']<=a['upper']))
        self.assertEqual(fit['parameters'][4],-4.214)
    def test_operator_reconstruction(self):
        from methanol_release.operators import verify_operators
        report=verify_operators();self.assertEqual(report['states'],204);self.assertLess(report['max_operator_error_khz'],1e-10)
    def test_independent_physical_hamiltonians(self):
        from methanol_zeeman.frozen_model import hyperfine
        from methanol_zeeman.spin import projection
        from methanol_hyperfine.spin import projection as native
        from methanol_zeeman.catalog import _inputs
        p=load_problem();a=np.load(DATA/'endpoint.npz',allow_pickle=False);h=p.d.hamiltonians(p.unpack(a['x']));waves=_inputs()['available_wavefunctions']
        for s,mat in h.items():
            overlap=projection(s.j,s.species)[0].conj().T@native(s.j,s.species)[0]
            other=hyperfine(s.key,waves[s.key]);self.assertLess(np.max(abs(overlap.conj().T@other@overlap-mat)),1e-9,s.key)
    def test_conditional_draws_at_all_noise_scales(self):
        from methanol_release.bootstrap import replay
        self.assertTrue(replay()['passed'])
    def test_magnetic_components_and_printed_cells(self):
        from methanol_zeeman.catalog import compute_transition
        source=ROOT/'supplement/zeeman/all_maser_components.csv';rows=list(csv.DictReader(source.open()))
        self.assertEqual(len(rows),750)
        for family in dict.fromkeys(r['family'] for r in rows):
            fresh=compute_transition(family,True);old=[r for r in rows if r['family']==family]
            self.assertEqual(len(fresh),len(old))
            for a,b in zip(fresh,old):
                for key in ['delta_nu_kHz','z_Hz_per_mG','g1_upper','g2_lower','statistical_line_fraction']:
                    self.assertAlmostEqual(a[key],float(b[key]),places=10)
        from methanol_release.tables import validate_zeeman_tables
        self.assertEqual(validate_zeeman_tables()['checked_rows'],224)
    def test_response_scenarios(self):
        from methanol_line_shift.estimate import load_components,estimate_shift
        rows=list(csv.DictReader((ROOT/'src/methanol_line_shift/data/profile_results.csv').open()))
        self.assertEqual(len(rows),232)
        for r in rows:
            f,s=load_components(r['transition']);fresh=estimate_shift(f,s,float(r['component_absorption_FWHM_kHz']),r['response'])
            self.assertLess(abs(fresh['shift_Hz']-float(r['fitted_center_Hz'])),1e-7)
            self.assertLess(abs(fresh['fitted_absorption_FWHM_kHz']-float(r['fitted_absorption_FWHM_kHz'])),1e-9)
    def test_experiment_family_regeneration(self):
        from methanol_release.experiment_families import calculate
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/'families';report=calculate(out)
            self.assertEqual(report['transitions'],128)
            self.assertLess(report['maximum_grid_position_error_kHz'],1e-4)
            self.assertEqual((out/'families.csv').read_bytes(),(ROOT/'supplement/experiment-families/families.csv').read_bytes())
    def test_reject_output_inside_source_and_existing(self):
        from methanol_zeeman.output_paths import validate_output_directory
        with self.assertRaises((ValueError,FileExistsError)):validate_output_directory(ROOT/'src/new-output',('a.csv',))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'existing';p.mkdir();(p/'a.csv').write_text('keep')
            with self.assertRaises((ValueError,FileExistsError)):validate_output_directory(p,('a.csv',))
            self.assertEqual((p/'a.csv').read_text(),'keep')
if __name__=='__main__':unittest.main()
