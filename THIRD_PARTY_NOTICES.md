# Scientific attribution and third-party notices

Project-authored software is licensed under EUPL 1.2. BELGI/GLOBALFIT,
the Xu constants, laboratory data and other third-party components retain
their own terms and the credits below. See [license scope](docs/LICENSING.md).

| Component or input | Attribution | Source |
|---|---|---|
| Torsion–rotation model and frozen constants | Li-Hong Xu, J. Fisher, R. M. Lees, H. Y. Shi, J. T. Hougen, J. C. Pearson, B. J. Drouin, G. A. Blake, R. Braakman | [Xu et al. 2008](https://doi.org/10.1016/j.jms.2008.03.017) |
| BELGI-Cs torsion–rotation program and method | J. T. Hougen, I. Kleiner, M. Godefroid | [Hougen, Kleiner and Godefroid 1994, J. Mol. Spectrosc. 163, 559–586](https://doi.org/10.1006/jmsp.1994.1047) |
| General methanol hyperfine theory and calculated tensors | Boy Lankhaar, Gerrit C. Groenenboom, Ad van der Avoird | [2016 paper](https://doi.org/10.1063/1.4972004), [2018 correction](https://doi.org/10.1063/1.5030801) |
| Magnetic tensors and Zeeman framework | Boy Lankhaar, Wouter Vlemmings, Gabriele Surcis, Huib Jan van Langevelde, Gerrit C. Groenenboom, Ad van der Avoird | [2018 paper](https://doi.org/10.1038/s41550-017-0341-8) |
| Measured rotational diagonal magnetic inputs | Kojiro Takagi, Shozo Tsunekawa, Kaori Kobayashi, Tomoya Hirota, Fusakazu Matsushima | [Takagi et al. 2021](https://doi.org/10.1016/j.jms.2021.111420) |
| Ground-state experimental profiles and comparison constants | L. H. Coudert, C. Gutlé, T. R. Huet, J.-U. Grabow, S. A. Levshakov | [Coudert et al. 2015](https://doi.org/10.1063/1.4926942) |
| Moderate/high-J splittings and triplet check | S. P. Belov and coauthors | [Belov et al. 2016](https://doi.org/10.1063/1.4954941) |
| Excited torsional-state profiles and tabular features | Li-Hong Xu and coauthors | [Xu et al. 2019](https://doi.org/10.1016/j.jms.2018.12.003) |
| Reproduced historical microwave profiles | J. E. M. Heuvel, A. Dymanus | [1973 paper I](https://doi.org/10.1016/0022-2852(73)90159-8), [paper II](https://doi.org/10.1016/0022-2852(73)90088-x) |
| Wigner numerical backend | C. Forssen, H. T. Johansson; Python source copyright Christian Forssen 2015 | [pywigxjpf 1.13.3](https://pypi.org/project/pywigxjpf/1.13.3/); Python header LGPL-3.0-or-later; package metadata GPLv3 |

The BELGI reference is “Selection Rules and Intensity Calculations for a Cs
Asymmetric Top Molecule Containing a Methyl Group Internal Rotor.” The methanol
Hamiltonian and parameter set are from Xu et al., “Torsion–rotation global
analysis of the first three torsional states (νt = 0, 1, 2) and terahertz
database for methanol.” Cite both for the corresponding program and inputs.
The retained BELGI/GLOBALFIT source lineage also credits Jonathan Fisher,
Ortigoso and Li-Hong Xu for its extensions; original source credits are retained.

The Wigner adapter is project-authored Python code wrapping the numerical
library. L26 uses
`pywigxjpf==1.13.3` for reproducible arithmetic; the magnetic contraction module
uses SymPy. The EUPL grant applies to the project-authored adapter, not to pywigxjpf
or the other dependencies. Attribution and terms for scientific inputs and
third-party components remain specific to those sources.

## Installed numerical libraries

| Library | Upstream license declaration | Project |
|---|---|---|
| NumPy | BSD | https://numpy.org |
| SciPy | BSD | https://scipy.org |
| SymPy | BSD | https://www.sympy.org |
| Matplotlib | Matplotlib license based on the PSF license | https://matplotlib.org |
| pywigxjpf 1.13.3 | GPLv3 package metadata; LGPL-3.0-or-later Python-wrapper header | https://pypi.org/project/pywigxjpf/1.13.3/ |

Consult each installed distribution for its complete license text and bundled
component notices. These external library distributions are not vendored in
this repository.
