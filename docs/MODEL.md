# Lankhaar (2026) model and reproducibility

L26 denotes the Lankhaar (2026) hyperfine model throughout this package.

The torsional momentum operator and wavefunctions use the native BELGI alpha
convention. Spatial tensors use exp(+i n alpha), with (a,b,c)=(z,x,y) and the
body-fixed commutator [Jb,Jc]=-i Ja. The DD source uses D_pi(i)pi(j)(-alpha),
where pi exchanges methyl sites 2 and 3. The ac, ca, bc and cb sine series of M
are reversed before fitting. Total-d source coefficients keep their declared
convention. The ordered M enters d-prime = d - M lambda; the displayed tensor
is its symmetric part. Lambda and rho are retained with the numerical inputs.

Fourier order is at most two. The 18 M series corrections lie in the affine
family `offset + M_basis @ amplitudes`; the nonzero offset is part of the
model. Four total-d constants and two DD scales complete the 13 molecular
fitting directions. The stored coordinate vector has 138 entries; its independent
representation has 127, including 114 response parameters. The fixed total-d
harmonic is -4.214 kHz. All bounds and the affine basis are packaged.

The objective partitions are Belov rows 0:55, retained Xu contrasts 55:99,
direct second-harmonic profiles 99:6432, and low-frequency profiles 6432:10142.
Low-frequency residuals are minimized subject to selective Q <= 931.127992775.
At L26, Q_selective = 931.1279927791534 and Q_low = 1590.0982142105986;
the 4.2e-9 guard excess is within the declared 1e-5 feasibility tolerance.
The total Q is 2521.2262069897515. These are half-sums of squared residuals,
with the retained data weights, masks, thermal background and crossover model.
Low-profile emission amplitudes use least squares followed by nonnegative
clipping without a second solve. This response rule is retained explicitly.

`methanol-release operators` reconstructs all 204 coefficient-operator sets
from theory tensors and eigenfunctions computed by bundled BELGI Fortran.
`fit --replay --derivatives` rebuilds the native operators and re-evaluates every
residual and derivative against the saved endpoint. Residuals must agree within
1e-10. Operator reconstruction changes rounding at about 1e-14 kHz; the profile
Jacobian uses finite differences, which amplify this rounding. Native derivatives
are compared within 1e-7 times max(1,abs(reference)). The separate
`fit --replay --derivatives --reference` check retains its original absolute
1e-10 derivative gate and exactly replays the saved operators.
`fit --run` performs a new constrained nonlinear optimization in the same
basis and bounds, starting at L26, into a separate directory. A bounded
restart does not overwrite or replace the packaged endpoint.

Uncertainties are conditional local-quadratic resampling ranges. Both stages
use the same seven-mode basis and shared selective noise. The zero-noise
selective minimum is 926.4872858631884; the fixed allowance is
4.640706911811662. The zero-noise surrogate displacement is removed before
propagation. The 768 archived draws comprise 512 nominal draws with a 1-kHz
table-noise assumption, 128 at 0.5 kHz and 128 at 2 kHz (PCG64, seed 20261019).
Printed uncertainties are the nominal 16th–84th percentiles. They condition on
the tensor family, basis, assignments, weights and noise assumptions; they are
not calibrated experimental confidence intervals or full nonlinear bootstrap
intervals. `bootstrap --replay` checks archived draws using their recorded
linearization. `bootstrap --run --output ...` recomputes the linearization from
native BELGI waves before generating new draws; optimizer roundoff need not give
byte-identical displacements. Neither command changes the published L26.

Magnetic inputs retain the measured rotational diagonal terms and computed
remaining tensor elements, with the declared single-c magnetic conversion.
Hyperfine eigenstates, component strengths and Zeeman factors are calculated
with L26 and native BELGI waves. Nominal line frequencies and dipole moments retain their stated
literature sources. Line-center shifts use full-real-axis least squares with
free amplitude, center and width, uniform frequency weights and zero baseline.
Statistical component populations are assumed; maser excitation is not modeled.

Table VIII uses a declared feature-selection rule: the finite-modulation 2f
phase has a positive isolated-line center; peaks must exceed 10% of the profile
maximum in height and 5% in prominence. Lorentz HWHM is 5 kHz and modulation
amplitude is 3 kHz, with statistical strengths and no crossovers. The counts
are computational descriptors, not calibrated experimental resolution limits.
All 128 transitions were checked at 0.02 and 0.01 kHz sampling, with refined
positions agreeing within 3.44e-7 kHz. Per-transition results also retain looser
and stricter thresholds. Reproduce with `python -m methanol_release.experiment_families
--output /absolute/new-family-results`.
