"""Regenerate response components from native BELGI waves and L26 tensors."""

import csv
import hashlib
import json

import numpy as np

from methanol_zeeman.catalog import DATA as MAGNETIC_DATA, _inputs
from methanol_zeeman.frozen_model import hyperfine, parameters
from methanol_zeeman.spin import transition_components
from methanol_zeeman.output_paths import validate_output_directory


def rebuild_response_catalog(output_dir, reconstruct_hyperfine=True):
    """Write all 136 response-transition patterns without fitting the molecule.

    Always recontracts physical tensors and native BELGI wavefunctions. The
    legacy reconstruct_hyperfine flag is accepted for API compatibility.
    """
    destination = validate_output_directory(output_dir,
        ("components.csv", "response_centroids.csv", "response_catalog_provenance.json"))
    inputs = _inputs()
    waves = inputs["available_wavefunctions"]
    matrices = {}
    source = MAGNETIC_DATA / "inputs/response_transition_list.json"
    transitions = {(r["upper"], r["lower"]): r["transition"] for r in json.loads(source.read_text())}
    rows, summaries = [], []
    for (upper, lower), name in transitions.items():
        for state in (upper, lower):
            if state not in matrices:
                matrices[state] = hyperfine(state, waves[state])
        ju, jl = int(upper.split("_")[2]), int(lower.split("_")[2])
        components = transition_components(ju, jl, upper[0], matrices[upper], matrices[lower], True)
        for offset, strength, fu, fl, pu, pl, iu, il in components:
            rows.append(dict(fit="L26", transition=name, upper=upper, lower=lower,
                offset_kHz=float(offset), strength=float(strength), F_upper=int(fu), F_lower=int(fl),
                PI_upper=int(pu), PI_lower=int(pl), upper_index=int(iu), lower_index=int(il)))
        summaries.append(dict(fit="L26", transition=name, upper=upper, lower=lower,
            components=len(components), centroid_Hz=float(np.average(components[:, 0], weights=components[:, 1]) * 1000),
            strength_sum=float(components[:, 1].sum()), span_kHz=float(np.ptp(components[:, 0]))))
    destination.mkdir(parents=True, exist_ok=True)
    for name, records in (("components.csv", rows), ("response_centroids.csv", summaries)):
        with (destination / name).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    report = dict(endpoint_sha256=parameters["source_sha256"], fit_performed=False,
        wavefunction_provider='native BELGI', hyperfine_recontracted=True, transitions=len(summaries), component_rows=len(rows),
        max_abs_centroid_Hz=max(abs(row["centroid_Hz"]) for row in summaries),
        source_transition_list_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    with (destination / "response_catalog_provenance.json").open("x") as handle:
        handle.write(json.dumps(report, indent=2) + "\n")
    return report
