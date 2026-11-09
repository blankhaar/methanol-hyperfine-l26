"""Validated public interface to the original full-real-axis line estimator.

Offsets and widths are in kHz. The returned shift is the fitted center minus
the statistical-strength centroid, in Hz. This is a specified measurement
operator, not a change to the underlying component frequencies or a maser model.
"""

import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np

from .profile_estimators import SHAPES, broad_limit, fitted_center

DATA = Path(__file__).resolve().parent / "data"
TARGETS = {"12.2GHz": "C2015-p06-1", "96.741GHz": "benchmark-96.741GHz"}


@lru_cache(maxsize=None)
def _native_components(transition):
    from methanol_zeeman.catalog import DATA as magnetic_data, _inputs
    from methanol_zeeman.frozen_model import hyperfine
    from methanol_zeeman.spin import transition_components
    records = json.loads((magnetic_data/'inputs/response_transition_list.json').read_text())
    row = next((r for r in records if r['transition']==transition), None)
    if row is None:
        raise ValueError(f'Unknown response transition {transition!r}')
    upper, lower = row['upper'], row['lower']
    waves = _inputs()['available_wavefunctions']
    return transition_components(int(upper.split('_')[2]), int(lower.split('_')[2]), upper[0],
        hyperfine(upper,waves[upper]), hyperfine(lower,waves[lower]), True)[:,:2]


def load_components(transition="12.2GHz", csv_path=None):
    """Calculate native BELGI/L26 components, or read an explicit user CSV.

    A supplied CSV must contain offset_kHz and strength; when it has a
    transition column, select the requested transition. No experimental file
    path or hidden workspace state is required.
    """
    name = TARGETS.get(transition, transition)
    if csv_path is None:
        components = _native_components(name)
        return components[:,0].copy(), components[:,1].copy()
    path = Path(csv_path)
    with path.open(newline="") as handle:
        records = list(csv.DictReader(handle))
    if records and "transition" in records[0]:
        records = [row for row in records if row["transition"] == name]
    if not records:
        raise ValueError(f"No components for transition {transition!r} in {path.name}")
    return (np.array([float(row["offset_kHz"]) for row in records]),
            np.array([float(row["strength"]) for row in records]))


def estimate_shift(offset_kHz, strength, width_kHz, response="Gaussian"):
    """Fit amplitude, center and width of one symmetric response to a blend.

    ``width_kHz`` is the component *absorption* FWHM for every response,
    including derivative responses. The minimization uses full-real-axis,
    uniformly weighted L2 with known zero baseline. The retained stationary
    search is designed for unresolved blends; a resolved pattern may have
    multiple competing optima. No global-optimum guarantee is made there.
    """
    offsets, weights = np.asarray(offset_kHz, dtype=float), np.asarray(strength, dtype=float)
    if offsets.ndim != 1 or weights.ndim != 1 or offsets.shape != weights.shape or offsets.size == 0:
        raise ValueError("Offsets and strengths must be nonempty, equal-length one-dimensional arrays.")
    if not np.isfinite(offsets).all() or not np.isfinite(weights).all():
        raise ValueError("Offsets and strengths must be finite.")
    if (weights < 0).any() or weights.sum() <= 0:
        raise ValueError("Strengths must be nonnegative with positive sum.")
    if not np.isfinite(width_kHz) or width_kHz <= 0:
        raise ValueError("Component absorption FWHM must be finite and positive.")
    if response not in SHAPES:
        raise ValueError(f"Choose response from {SHAPES}.")
    weights = weights / weights.sum()
    centroid = float(weights @ offsets)
    result = fitted_center(offsets, weights, float(width_kHz), response)
    return dict(response=response, component_absorption_FWHM_kHz=float(width_kHz),
        statistical_centroid_kHz=centroid, fitted_center_kHz=float(result["center"]),
        shift_Hz=float(1000 * (result["center"] - centroid)),
        fitted_absorption_FWHM_kHz=float(result["fitted_absorption_FWHM"]),
        amplitude=float(result["amplitude"]), stationarity_norm=result["stationarity_norm"],
        broad_limit_Hz=float(1000 * broad_limit(offsets, weights, width_kHz, response)))
