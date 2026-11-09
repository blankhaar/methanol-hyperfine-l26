# Computing the torsion–rotation eigenfunctions

The default provider compiles and runs the bundled BELGI-Cs / methanol GLOBALFIT
Fortran code. It supplies the eigenfunctions to the hyperfine operator builder,
L26 likelihood and optimizer, magnetic catalogues, line-shift calculations and
Table VIII. No separate research checkout or precomputed eigenfunction archive
is required for these calculations.

## Source and citation

The 58 upstream files are retained byte for byte under
`src/methanol_hyperfine/data/belgi/source/`, including their original credits.
The adjacent `PROVENANCE.json` and `source/SHA256SUMS` identify the exact source
snapshot. The program is credited to J. T. Hougen, I. Kleiner and M. Godefroid,
[J. Mol. Spectrosc. 163, 559–586 (1994)](https://doi.org/10.1006/jmsp.1994.1047).
The methanol Hamiltonian and its constants are from
[Xu et al., J. Mol. Spectrosc. 251, 305–313 (2008)](https://doi.org/10.1016/j.jms.2008.03.017).
Cite both works. Original Fisher, Ortigoso and Xu extension credits are retained.

The 119 named assignments in `inputs/xu2008-final-constants.txt` are mapped to
the complete 205-entry BELGI CST vector; unassigned entries are zero. The parser
checks the exact constants checksum, duplicate/unknown names and finite values.
The torsion–rotation constants are fixed. They are separate from the fitted
hyperfine coupling constants of L26.

## Build and execution

Install the Python package and GNU Fortran. On macOS the build links Apple
Accelerate; on Linux install LAPACK and BLAS development libraries. This native
adapter supports macOS and Linux. Set `METHANOL_FORTRAN_COMPILER` to select an
explicit GNU Fortran executable if it is not named `gfortran` on PATH.

```sh
python -m pip install -e .
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
methanol-release belgi --force --output /absolute/fresh-belgi-run
methanol-release fit --replay --derivatives
python reproduce.py --output /absolute/fresh-reproduction
```

The first command computes the 204 states in the fit. Individual states can be
requested repeatedly with `--state`, for example `--state A_0_5_1_+` or
`--state E_1_6_-2_E`. The driver supports A/E states with J <= 46, vt = 0, 1, 2,
and 4 <= ktronc <= 20; L26 uses ktronc = 10. These bounds describe the interface,
not a claim of validated spectroscopic accuracy over every possible state.
`--all-reference-states --ktronc 10` computes the full reference state list.

The driver uses BELGI's ENCAL diagonalization, state tracking and expansion of
the contracted torsional vectors into the original Fourier basis. Energies are
in cm^-1; coefficient rows are K = -J,...,J and columns n = -ktronc,...,ktronc
with m = 3n + sigma. It uses native alpha and p_alpha, state-tracking label 1,
and no post hoc coefficient renormalization. The existing ENCAL behavior that
zeros the ninth contracted torsional vector is retained and recorded. The
stored reference waves use the same convention and limitation.

One existing array-bound compatibility fix is applied to a temporary build
copy of SETUP: `max = 2*nn2+1` becomes `max = kdmx`, preventing a 95-element clear
of 93-element arrays. Original vendored files remain unchanged. The project
driver's request capacity is 13254, sufficient for all supported state labels.
The build receipt records compiler flags, source hashes and this exact patch.

## Generated caches and reference checks

Builds and computed waves are stored outside the package, by default under
`~/.cache/methanol-hyperfine/belgi`. Set `METHANOL_BELGI_CACHE` to a different
directory or use `belgi --cache-dir`. Cache identities include source, constants,
driver, interface, compiler binary/version and platform hashes. Result receipts
include requested state identities, truncation, energies, norms and array hashes.
`--force` recalculates waves. Missing compilers, changed scientific inputs and
damaged generated caches fail explicitly; stored reference waves are never an
automatic fallback.

Reference waves and operator matrices remain available to check regeneration.
`export_wavefunctions(..., backend='reference')` and
`fit --replay --derivatives --reference` select those checks explicitly. The
normal `fit --replay` uses newly constructed native operators. Their rounding
can differ at floating-point precision; see [MODEL.md](MODEL.md) for separate
residual and finite-difference Jacobian tolerances. Neither route refits or
changes the packaged endpoint.

`bootstrap --replay`, `shift --surrogate`, and the figure recipes reproduce
recorded derived products. New bootstrap runs and ordinary line-shift estimates
use native BELGI. Rebuilding response products also recalculates the components
before constructing the bounded surrogates. Published tables, figures and
conditional uncertainty ranges use the L26 values.
