"""Fast approximations to the frozen statistical-strength line-pull curves.

Input W is the individual component absorption FWHM in kHz, even for
derivative responses. Output is fitted center minus statistical centroid, Hz.
Valid only on 30 <= W <= 2000 kHz and for the archived response/fit assumptions.
The approximation error is not an uncertainty on the molecular prediction.
"""
from functools import lru_cache
import json
from pathlib import Path

import numpy as np

COEFFICIENTS = Path(__file__).resolve().parent / "data" / "line_pull_surrogate_coefficients.json"


@lru_cache(maxsize=1)
def coefficient_catalog():
    return json.loads(COEFFICIENTS.read_text())


def line_pull_Hz(width_kHz, target="12.2GHz", response="Gaussian"):
    """Evaluate a rational approximation; reject unsupported extrapolation.

    Targets: '12.2GHz', '96.741GHz'. Responses: 'Gaussian', 'Lorentz',
    'dispersion_1f', 'absorption_2f'. Accepts a scalar or an array.
    """
    width = np.asarray(width_kHz, dtype=float)
    if np.any(~np.isfinite(width)) or np.any((width < 30.) | (width > 2000.)):
        raise ValueError("Component absorption FWHM must lie in [30, 2000] kHz.")
    key = f"{target}/{response}"
    models = coefficient_catalog()["models"]
    if key not in models:
        raise ValueError(f"Unknown target/response {key}; choose from {list(models)}")
    model = models[key]
    x = (30. / width)**2
    numerator = np.polynomial.polynomial.polyval(x, model["numerator"])
    denominator = np.polynomial.polynomial.polyval(x, model["denominator"])
    result = model["A_Hz_kHz2"] / width**2 * numerator / denominator
    return float(result) if width.ndim == 0 else result
