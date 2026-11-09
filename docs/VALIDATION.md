# Numerical validation

The tests compare L26 calculations against the packaged numerical references.
The reference arrays support regression checks; the default numerical provider
computes eigenfunctions using the bundled BELGI source.

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
