# Numerical validation

The tests compare L26 calculations against the packaged numerical references.
The reference arrays support regression checks; the default numerical provider
computes eigenfunctions using the bundled BELGI source.

## Complete paper reconstruction

`scripts/reproduce_paper.py` connects the numerical calculations to the table
and figure renderers. See the installation and two-renderer commands in the
[README](../README.md). Its default run uses a new native BELGI cache inside
the chosen output and recalculates all 768 conditional draws. The eight main
tables and 27 supplementary magnetic tables are written as publication-ready
TeX fragments, along with a machine-readable coefficient/interval table.
All 34 profile predictions, 57 Belov splittings, the independent triplet
profile and 1,464 line-pull ordinates are recalculated before rendering.

The retained L26 endpoint, response parameters, molecular input tensors,
observations, assignments and uncertainty prescription define this calculation.
The figure and table reconstruction does not run a new molecular optimization,
repeat the electronic-structure calculation or redigitize published graphs.
Literature comparisons, captions and display precision are retained inputs.
Reference plotting arrays are read only after fresh predictions have been made.

The run requires exact agreement with the printed TeX tables, profile and
Belov differences below 1e-9 in their respective signal/kHz units, and line-pull
differences below 1e-6 Hz. PDF hashes are reported separately: numerical agreement
is meaningful across platforms even when renderer or font versions change PDF
bytes. `validation.json` and the per-step logs record each check.

## Regression checks

```sh
python scripts/run_tests.py
methanol-release fit --replay --derivatives
methanol-release operators
methanol-release bootstrap --replay
python scripts/verify_manifest.py --source-tree
```

The validation records under `release-evidence/` cover the following checks:

- All 4,148 reference eigenfunctions across three truncations regenerate on
  macOS with GNU Fortran and Apple Accelerate.
- The 10,142 fitted residuals agree to about 1.3e-13. Native finite-difference
  derivatives satisfy the declared 1e-7 scaled tolerance; the exact stored-operator
  replay has a separate 1e-10 absolute gate.
- Independent Hamiltonians cover 204 states. Magnetic tables contain 750
  components for 27 transition families; all 224 printed rows are checked.
- Conditional resampling has 512 nominal draws and two 128-draw noise controls.
  These are local approximation checks, not uncertainty calibration.
- Candidate experiments cover 128 transitions. The feature-selection rule and
  sampling/threshold checks are specified in [MODEL.md](MODEL.md).
- Installation checks exercise fitting, magnetic and line-shift calculations
  with the packaged reference wave/operator archives removed.

Linux linking is supported by the build adapter, but Linux numerical execution
has not been validated in these records. The compiler, platform and reference
identities attached to individual receipts define their tested scope. Validation
of a numerical implementation does not establish physical identifiability,
calibrated uncertainties or maser excitation predictions.
