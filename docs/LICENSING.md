# License and third-party components

Copyright (c) 2026 Boy Lankhaar. The project-authored software is licensed under
**EUPL 1.2**; the complete text is in [LICENSE](../LICENSE). The SPDX identifier is
`EUPL-1.2`. The [European Commission's license page](https://interoperable-europe.ec.europa.eu/collection/eupl/eupl-text-eupl-12)
provides the official versions.

The grant covers the Python code under `src/`, the numerical, reproduction and
figure-generation scripts, the tests, and the project prediction driver
`src/methanol_hyperfine/data/belgi/fortran/belgi_wavefunctions.f`.
The scientific input data and publication assets are outside this software
license grant. The exceptions below retain their own terms and attribution.

## BELGI and molecular constants

The 58 original BELGI/GLOBALFIT Fortran files in
`src/methanol_hyperfine/data/belgi/source/` are third-party software, retained
byte for byte with their original author notices. The Xu constants in the
adjacent `inputs/` directory are separately attributed scientific inputs.
Neither is relicensed under EUPL by this project.

BELGI-Cs is credited to J. T. Hougen, I. Kleiner and M. Godefroid,
[J. Mol. Spectrosc. 163, 559–586 (1994)](https://doi.org/10.1006/jmsp.1994.1047).
The methanol Hamiltonian and parameters are credited to
[Xu et al., J. Mol. Spectrosc. 251, 305–313 (2008)](https://doi.org/10.1016/j.jms.2008.03.017).
The source also retains the Fisher, Ortigoso and Xu extension credits. Cite
both program and molecular-model references when using these calculations.
Source hashes, the build compatibility patch and numerical conventions are
documented in [BELGI.md](BELGI.md).

## Experimental inputs and dependencies

Digitized laboratory spectra and literature tensor inputs retain their source
attribution; the software license does not relicense the underlying articles,
supplements or third-party data. [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md)
identifies the publications and numerical libraries.

NumPy, SciPy, SymPy, Matplotlib and pywigxjpf are installed as dependencies and
retain their own licenses. Their distributions are not bundled here. In
particular, pywigxjpf 1.13.3 declares GPLv3 in its package metadata while its
Python-wrapper source header declares LGPL-3.0-or-later. Both upstream notices
are preserved as attribution; this package does not replace either declaration.

The explicit export inventory and checksums define this source distribution.
The export checks reject known excluded source documents and dependency binaries;
they do not replace the component-specific license notices above.
