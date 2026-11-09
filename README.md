# L26 methanol hyperfine and Zeeman spectroscopy

This package implements the **Lankhaar (2026) model**, abbreviated **L26**, for
methanol hyperfine spectra, weak-field Zeeman coefficients and line-center
shifts. It provides the molecular tensors, fitting routines, conditional
uncertainty ranges, machine-readable catalogues, and the inputs and recipes
for the paper's five figures and seven supplementary plates.

Torsion–rotation eigenfunctions are computed from the bundled BELGI Fortran
program and Xu molecular constants. Published Belov Fit A(2) is included as a
literature comparison in Figure 2.

## Install

Use Python 3.12 or later and GNU Fortran (`gfortran`). On macOS, the native build
links Apple Accelerate. Linux requires LAPACK and BLAS development libraries.
The exact tested numerical environment is in `requirements-numerical.txt`.

```sh
git clone https://github.com/blankhaar/methanol-hyperfine-l26.git
cd methanol-hyperfine-l26
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-numerical.txt
python -m pip install -e .
```

The first calculation compiles BELGI and computes the required eigenfunctions.
Subsequent calls may reuse locally generated results with matching source,
constants and compiler fingerprints. Set `METHANOL_FORTRAN_COMPILER` to choose
another GNU Fortran executable; [BELGI.md](docs/BELGI.md) describes the build
and cache controls.

## Calculate and reproduce

Use one numerical thread for reproducible replay:

```sh
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
methanol-release fit --replay --derivatives
methanol-release operators
methanol-release bootstrap --replay
methanol-release zeeman --help
methanol-release shift --help
```

To recompute BELGI eigenfunctions or run a constrained nonlinear fit:

```sh
methanol-release belgi --force --output /absolute/new-belgi-waves
methanol-release fit --run --output /absolute/new-fit-run
```

L26 has four fitted total-d constants, seven affine spin-rotation directions
and two dipole-dipole scales, together with experimental-response parameters.
The complete nonlinear forward model evaluates 10,142 residuals.
[MODEL.md](docs/MODEL.md) specifies the conventions, bounds, objective and
conditional uncertainty method.

The reported uncertainty ranges condition on the tensor family, basis,
assignments, weights and noise assumptions. They use local-linear and
local-quadratic approximations and are not calibrated experimental confidence
intervals or a full nonlinear bootstrap.

## Figures, tables and catalogues

- `current-figures/`: twelve publication figures, numerical inputs and rendering scripts.
- `tables/`: numerical LaTeX table fragments, including literature comparisons.
- `supplement/`: coupling coefficients, signed modes, conditional uncertainty ranges,
  magnetic components, line shifts and candidate experimental transitions.

After installation, the following command rebuilds the eigenfunctions,
spectra, all 768 conditional uncertainty draws, eight main tables, 27
supplementary tables, and all twelve figures. It recomputes the plotted model
curves; saved plotting arrays are used only for comparison. The public L26
endpoint and experimental-response parameters are inputs, so this reproduces
the published model without performing a new molecular fit.

The two rendering environments preserve the paper's original appearance:
main figures use Matplotlib 3.11.2 and supplementary plates use 3.10.7.
With the numerical environment still activated:

```sh
python -m venv ../l26-main-figures-venv
../l26-main-figures-venv/bin/python -m pip install numpy==2.3.5 matplotlib==3.11.2
python scripts/reproduce_paper.py --output /absolute/new-paper-run \
  --main-python "$PWD/../l26-main-figures-venv/bin/python" \
  --atlas-python "$PWD/.venv/bin/python"
```

Results are in `figures/manuscript/figures/`, `figures/supplement/figures/`
and `tables/` below the chosen output. `validation.json` reports the numerical,
table and figure comparisons. The table builder retains publication layouts,
display precision and literature values while recalculating every L26/theory
cell and conditional interval. No private manuscript or local research files
are needed. A full run takes several minutes, depending on hardware.

For a quicker redraw using the packaged plotting arrays, use
`current-figures/reproduce.py` with the same renderer options. `reproduce.py`
alone runs the numerical regression/catalogue checks, including three reference
uncertainty draws; it does not rebuild the complete paper products.

Output directories must be fresh and outside the source tree.
[VALIDATION.md](docs/VALIDATION.md) describes the numerical checks and their limits.

## Citation and license

Please cite the accompanying Lankhaar (2026) paper and the relevant input sources
listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
The torsion–rotation calculation uses **BELGI-Cs**, credited to J. T. Hougen,
I. Kleiner and M. Godefroid
([1994, J. Mol. Spectrosc. 163, 559–586](https://doi.org/10.1006/jmsp.1994.1047)),
and the methanol Hamiltonian and constants of
[Xu et al. (2008, J. Mol. Spectrosc. 251, 305–313)](https://doi.org/10.1016/j.jms.2008.03.017).

Project-authored software is licensed under **EUPL 1.2**; see [LICENSE](LICENSE).
BELGI/GLOBALFIT, the Xu constants, experimental inputs and numerical dependencies
retain their own terms and attribution. The precise scope is in
[LICENSING.md](docs/LICENSING.md).
