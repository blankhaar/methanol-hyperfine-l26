"""Compute magnetic components from native BELGI waves and fixed L26 tensors.

The default rebuilds hyperfine eigenstates and transition strengths. Packaged
catalogues remain comparison references, not a substitute eigenfunction provider.
"""

import csv
import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.constants import c, epsilon_0, h

from .contractions import Contractions, fourier_tensor
from .frozen_model import hyperfine, parameters, wave
from methanol_hyperfine.belgi import State, export_wavefunctions
from .magnetic_angular import MU_N_HZ_PER_MG, PROTON_G, weak_level_g, weak_transition_g
from .spin import diagonalize, transition_components
from .wigner_adapter import wig3jj
from .output_paths import validate_output_directory

DATA = Path(__file__).resolve().parent / "data"


def _read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


@lru_cache(maxsize=1)
def _inputs():
    """Build native waves from transition identities, without opening saved waves."""
    result = {}
    result["transitions"] = json.loads((DATA / "inputs/transitions.json").read_text())["transitions"]
    responses = json.loads((DATA / "inputs/response_transition_list.json").read_text())
    states = {State.from_key(row[end]) for row in [*result['transitions'], *responses]
              for end in ('upper', 'lower')}
    result['available_wavefunctions'] = {s.key: w.coefficients for s, w in export_wavefunctions(states).items()}
    result["magnetic"] = json.loads((DATA / "inputs/magnetic_tensor_input.json").read_text())
    return result


def contract_magnetic(state, coefficients=None, gx=None):
    """Return rotation, effective torsion and scalar-reduction checks for a level.

    ``coefficients`` maps aa,...,cc to Fourier coefficients at n=0,3,6.
    The native alpha convention requires the single-c sine sign conversion
    exactly once below. gx=0.425 is already effective: do not apply the separate
    g_alpha to gx conversion a second time.
    """
    data = _inputs()
    definition = data["magnetic"]
    coefficients = definition["magnetic_inputs"]["g_Fourier"] if coefficients is None else coefficients
    gx = definition["adoption"]["effective_torsion_gx"] if gx is None else gx
    coeff = data["available_wavefunctions"][state]
    con = Contractions(wave(state, coeff))
    if con.jj == 0:
        raise ValueError("The scalar rotational contraction requires J > 0.")
    ns = 1 if state.startswith("A_") else 2
    rotation = np.zeros((ns, ns), complex)
    for index, harmonic in enumerate((0, 3, 6)):
        cosine, sine = np.zeros((3, 3)), np.zeros((3, 3))
        for pair, values in coefficients.items():
            a, b = map("abc".index, pair)
            (sine if (a == 2) != (b == 2) else cosine)[a, b] = values[index]
        rotation += con.tensor(fourier_tensor(cosine, -sine, harmonic)) / con.jj
    torsion = gx * con.spin_torsion(0, 0) / con.jj
    total = rotation + torsion
    scalar_error = float(np.max(abs(total - np.eye(ns) * np.trace(total) / ns)))
    if scalar_error >= 1e-12:
        raise ArithmeticError(f"Magnetic E-partner reduction is not scalar for {state}.")
    weights = [float(con.quadratic(a, a, 0)[0, 0].real / con.jj) for a in range(3)]
    return {
        "rotation": float(rotation[0, 0].real),
        "torsion": float(torsion[0, 0].real),
        "gTR": float(total[0, 0].real),
        "partner_scalar_error": scalar_error,
        "diagonal_weights_abc": weights,
    }


def _state_levels(state, rebuild=True, reconstruct_hyperfine=True):
    data = _inputs()
    j, species = int(state.split("_")[2]), state[0]
    waves = data["available_wavefunctions"]
    matrix = hyperfine(state, waves[state])
    energy, vectors, fs, parity = diagonalize(j, species, matrix)
    magnetic = contract_magnetic(state)
    g_values = weak_level_g(j, species, vectors, fs, magnetic["gTR"])
    a_j = weak_level_g(j, species, vectors, fs, 1., gp=0.)
    rows = []
    for i, value in enumerate(energy):
        siblings = sorted(np.flatnonzero(fs == fs[i]), key=lambda k: (energy[k], parity[k], k))
        label = str(int(fs[i])) + (chr(ord("a") + siblings.index(i)) if len(siblings) > 1 else "")
        rows.append(dict(state=state, index=i, F=int(fs[i]), PI=int(parity[i]), label=label,
            hyperfine_energy_kHz=float(value), A_J=float(a_j[i]), g=float(g_values[i]),
            gTR=magnetic["gTR"], g_rotation=magnetic["rotation"], g_torsion=magnetic["torsion"],
            partner_scalar_error=magnetic["partner_scalar_error"]))
    return rows, matrix


def _magnetic_columns(row, upper, lower):
    fu, fl = int(row["F_up"]), int(row["F_down"])
    geff = weak_transition_g(upper["g"], lower["g"], fu, fl)
    return dict(z_Hz_per_mG=float(-2 * MU_N_HZ_PER_MG * geff),
                g1_upper=upper["g"], g2_lower=lower["g"], g_eff=float(geff))


def compute_transition(family, reconstruct_hyperfine=False):
    """Rebuild a catalogued transition's hyperfine strengths and Zeeman response.

    Nominal frequencies and dipole moments determine Einstein A, not hyperfine
    offsets. H is always contracted from L26 tensors and native BELGI waves.
    The legacy reconstruct_hyperfine flag is accepted for API compatibility.
    """
    data = _inputs()
    targets = {t["label_GHz"]: t for t in data["transitions"]}
    if family not in targets:
        raise ValueError(f"Unknown family {family!r}; choose from {list(targets)}")
    target = targets[family]
    upper, lower = target["upper"], target["lower"]
    ju, jl, species = int(upper.split("_")[2]), int(lower.split("_")[2]), upper[0]
    upper_levels, hu = _state_levels(upper, True, reconstruct_hyperfine)
    lower_levels, hl = _state_levels(lower, True, reconstruct_hyperfine)
    components = transition_components(ju, jl, species, hu, hl, True)
    summary = next(r for r in _read_csv(DATA / "reference/transition_summary.csv") if r["family"] == family)
    frequency_hz = float(summary["nominal_frequency_MHz"]) * 1e6
    cu, cl = (data["available_wavefunctions"][key] for key in (upper, lower))
    amplitude = 0j
    # One E partner's reduced rotational amplitude. Partner multiplicity belongs
    # in the hyperfine strength sum, not in A per upper level.
    for q, dipole in {0: .899, 1: 1.44 / np.sqrt(2), -1: -1.44 / np.sqrt(2)}.items():
        for iu, ku in enumerate(range(-ju, ju + 1)):
            kl = ku - q
            if not -jl <= kl <= jl:
                continue
            il = kl + jl
            amplitude += (np.vdot(cu[iu], cl[il]) * (-1.) ** (ju - ku)
                * np.sqrt((2 * ju + 1) * (2 * jl + 1))
                * wig3jj(2 * ju, 2, 2 * jl, -2 * ku, 2 * q, 2 * kl) * dipole)
    smu2 = float(abs(amplitude) ** 2)
    a_rot = 16 * np.pi**3 * frequency_hz**3 / (3 * epsilon_0 * h * c**3) * (1e-21 / c)**2 * smu2 / (2 * ju + 1)
    total = float(components[:, 1].sum())
    centroid = float(np.average(components[:, 0], weights=components[:, 1]))
    rows = []
    for off, strength, fu, fl, pu, pl, iu, il in components:
        fu, fl, iu, il = map(int, (fu, fl, iu, il))
        u, l = upper_levels[iu], lower_levels[il]
        relative_a = (2 * ju + 1) * strength / (2 * fu + 1)
        row = dict(family=family, upper_state=upper, lower_state=lower,
            nominal_frequency_MHz=frequency_hz / 1e6, F_up=fu, F_down=fl,
            upper_label=u["label"], lower_label=l["label"], upper_PI=int(pu), lower_PI=int(pl),
            upper_index=iu, lower_index=il, delta_nu_kHz=float(off - centroid))
        row.update(_magnetic_columns(row, u, l))
        row.update(log10_A_s_inverse=float(np.log10(a_rot * relative_a)),
            A_s_inverse=float(a_rot * relative_a), reduced_strength=float(strength),
            statistical_line_fraction=float(strength / total), relative_A=float(relative_a),
            z_coefficient_gTR_upper=float(-2 * MU_N_HZ_PER_MG * weak_transition_g(u["A_J"], 0, fu, fl)),
            z_coefficient_gTR_lower=float(-2 * MU_N_HZ_PER_MG * weak_transition_g(0, l["A_J"], fu, fl)),
            z_coefficient_gp=float(-2 * MU_N_HZ_PER_MG * weak_transition_g(1-u["A_J"], 1-l["A_J"], fu, fl)))
        rows.append(row)
    return rows


def calculate_catalog(output_dir, families=None, rebuild_components=False, reconstruct_hyperfine=False):
    """Write components, levels and provenance; return the run summary.

    Hyperfine eigenstates and components are rebuilt from native BELGI waves
    and L26 tensors. Legacy rebuild flags are accepted for API compatibility.
    """
    output_dir = validate_output_directory(output_dir,
        ("components.csv", "hyperfine_levels.csv", "provenance.json"))
    rebuild_components = True
    reconstruct_hyperfine = True
    reference = _read_csv(DATA / "reference/all_maser_components.csv")
    available = list(dict.fromkeys(row["family"] for row in reference))
    selected = available if families is None else list(dict.fromkeys(families))
    if not selected or not set(selected) <= set(available):
        raise ValueError(f"Choose at least one family from {available}")
    reference = [row for row in reference if row["family"] in selected]
    states = sorted({row[key] for row in reference for key in ("upper_state", "lower_state")})
    levels = {state: _state_levels(state, rebuild_components, reconstruct_hyperfine)[0] for state in states}
    rows = [row for family in selected for row in compute_transition(family, reconstruct_hyperfine)]
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, records in (("components.csv", rows), ("hyperfine_levels.csv", [r for state in states for r in levels[state]])):
        with (output_dir / name).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    report = dict(fit="L26", torsion_rotation_provider="native BELGI", endpoint_sha256=parameters["source_sha256"],
        fit_performed=False, eigenstates_recomputed=rebuild_components,
        hyperfine_recontracted=reconstruct_hyperfine, families=selected,
        components=len(rows), states=len(states),
        input_sha256={str(path.relative_to(DATA)): hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in sorted((DATA / "inputs").iterdir()) if path.is_file()},
        magnetic_definition=_inputs()["magnetic"],
        nuclear_moment=dict(mu_N_Hz_per_mG=MU_N_HZ_PER_MG, proton_g=PROTON_G, source="scipy.constants"),
        convention="E/h=E0/h-(mu_N/h)gM B; z=(nu_qplus-nu_qminus)/B=-2(mu_N/h)g_eff; q=Mu-Ml",
        printed_manuscript_tables="L26 numerical tables")
    with (output_dir / "provenance.json").open("x") as handle:
        handle.write(json.dumps(report, indent=2) + "\n")
    return report
