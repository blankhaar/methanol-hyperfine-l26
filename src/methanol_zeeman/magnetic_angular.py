"""Independent angular Zeeman utilities for an isolated torsion-rotation level.

No magnetic tensor, tensor-frame mapping, or wavefunction provider is assumed.
Inputs are eigenvectors returned by methanol_zeeman.spin.diagonalize at MF=0,
and a separately established scalar or partner-space magnetic g matrix.
These functions do not read constants, run BELGI, or write reference products.
"""
from functools import lru_cache
import itertools
import numpy as np
from scipy.constants import physical_constants, h
from .spin import spin_basis, projection
from .wigner_adapter import wig3jj

MU_N_HZ_PER_MG = physical_constants['nuclear magneton'][0] / h * 1e-7
PROTON_G = physical_constants['proton g factor'][0]


@lru_cache(None)
def uncoupled_basis(j, mf):
    """Spin order matches production: methyl protons first, hydroxyl last."""
    s = np.array(list(itertools.product([-.5, .5], repeat=4)))
    s = s[abs(mf-s.sum(axis=1)) <= j]
    return s, mf-s.sum(axis=1)


@lru_cache(None)
def total_step(j, mf, direction):
    """F+ or F- in explicit spin/rotor basis, all angular momenta in hbar."""
    assert direction in (-1, 1)
    s, mj = uncoupled_basis(j, mf)
    t, _ = uncoupled_basis(j, mf+direction)
    index = {tuple(x): i for i, x in enumerate(t)}
    out = np.zeros((len(t), len(s)))
    for col, x in enumerate(s):
        if tuple(x) in index:
            out[index[tuple(x)], col] += np.sqrt(max(0., j*(j+1)-mj[col]*(mj[col]+direction)))
        for proton in range(4):
            if x[proton] == -.5*direction:
                y = x.copy()
                y[proton] += direction
                if tuple(y) in index:
                    out[index[tuple(y)], col] += 1.
    return out


def eigenvectors_at_mf(j, species, vectors, fs, mf):
    """Return uncoupled vectors and their original eigenstate indices at MF."""
    assert int(mf) == mf
    ids = np.flatnonzero(fs >= abs(mf))
    fs = np.asarray(fs)[ids]
    out = projection(j, species)[0] @ vectors[:, ids]
    ns = 1 if species == 'A' else 2
    direction = 1 if mf >= 0 else -1
    for m in range(0, mf, direction):
        step = np.kron(np.eye(ns), total_step(j, m, direction))
        out = (step @ out) / np.sqrt(fs*(fs+1)-m*(m+direction))[None, :]
    return out, ids


def weak_level_g(j, species, vectors, fs, gtr, gp=PROTON_G):
    """g_F in E(B)/h = E(0)/h - (mu_N/h) B M_F g_F.

    gtr is required to be the same scalar in both E partner sectors. A matrix
    cannot silently be substituted. F=0 has zero linear magnetic shift.
    """
    if np.ndim(gtr) != 0:
        raise ValueError('Verify scalar E-partner reduction before calling')
    _, _, _, pairs, *_ = spin_basis(j)
    ii = 3*np.eye(next(iter(pairs.values())).shape[0]) + 2*sum(pairs.values())
    ns = 1 if species == 'A' else 2
    t = projection(j, species)[0]
    projected = t.conj().T @ np.kron(np.eye(ns), ii) @ t
    i2 = np.einsum('ai,ab,bi->i', vectors.conj(), projected, vectors).real
    f2 = fs*(fs+1)
    a_j = np.zeros(len(fs))
    nonzero = f2 > 0
    a_j[nonzero] = (f2[nonzero] + j*(j+1) - i2[nonzero])/(2*f2[nonzero])
    out = np.zeros(len(fs))
    out[nonzero] = gtr*a_j[nonzero] + gp*(1-a_j[nonzero])
    return out


def magnetic_matrix_at_mf(j, species, vectors, fs, mf, gtr, gp=PROTON_G):
    """Dimensionless magnetic moment matrix G: HZ/h = -mu_N/h * B * G.

    Matrix is in the surviving MF-zero eigenstate labeling, lifted to MF.
    gtr may be a scalar or a Hermitian ns-by-ns matrix (ns=1 for A, 2 for E).
    A generic matrix is projected onto the physical permutation subspace;
    callers remain responsible for its actual molecular symmetry and origin.
    Does not supply the rotational-state coupling or diamagnetism required
    for high-field laboratory modeling.
    """
    u, ids = eigenvectors_at_mf(j, species, vectors, fs, mf)
    spins, mj = uncoupled_basis(j, mf)
    ns = 1 if species == 'A' else 2
    gtr = np.asarray(gtr)
    if gtr.ndim == 0:
        gtr = complex(gtr)*np.eye(ns)
    if gtr.shape != (ns, ns) or not np.allclose(gtr, gtr.conj().T, atol=1e-12, rtol=0):
        raise ValueError('gtr must be scalar or a Hermitian E-partner matrix')
    operator = np.kron(gtr, np.diag(mj))
    operator += gp*np.kron(np.eye(ns), np.diag(spins.sum(axis=1)))
    g = u.conj().T @ operator @ u
    return (g+g.conj().T)/2, ids


def weak_transition_g(gu, gl, fu, fl):
    """Strength-weighted sigma centroid g, not a finite-field profile peak."""
    return (gu+gl)/2 + (gu-gl)*(fu*(fu+1)-fl*(fl+1))/4


def dipole_at_mf(ju, jl, species, vu, fu, vl, fl, mu, ml, partner_dipole=None):
    """Normalized q=Mu-Ml electric dipole between lifted zero-field bases.

    Default partner matrix matches production A/E electric-dipole convention.
    A caller with explicitly integrated reduced rotational amplitudes can pass
    those as partner_dipole instead. The two must not be applied twice.
    Returns matrix and original upper/lower state indices.
    """
    q = mu-ml
    if q not in (-1, 0, 1):
        raise ValueError('Electric-dipole q must be -1, 0 or +1')
    u, iu = eigenvectors_at_mf(ju, species, vu, fu, mu)
    l, il = eigenvectors_at_mf(jl, species, vl, fl, ml)
    su, mju = uncoupled_basis(ju, mu)
    sl, mjl = uncoupled_basis(jl, ml)
    lookup = {tuple(s): i for i, s in enumerate(sl)}
    d = np.zeros((len(su), len(sl)))
    for row, spin in enumerate(su):
        if tuple(spin) in lookup:
            col = lookup[tuple(spin)]
            d[row, col] = (-1.)**(ju-int(mju[row]))*wig3jj(
                2*ju, 2, 2*jl, -2*int(mju[row]), 2*q, 2*int(mjl[col]))
    ns = 1 if species == 'A' else 2
    if partner_dipole is None:
        partner_dipole = np.diag([1.] if ns == 1 else [1., -1.])
    if np.shape(partner_dipole) != (ns, ns):
        raise ValueError('Incorrect torsion-rotation partner dipole shape')
    return u.conj().T @ np.kron(partner_dipole, d) @ l, iu, il


def finite_field_block(j, species, energies_khz, vectors, fs, mf, gtr, field_mg,
                       gp=PROTON_G):
    """Isolated-J linear-Zeeman Hamiltonian including hyperfine-state mixing.

    Returns energies, zero-field-coordinate eigenvectors and original indices.
    Ordinary Hermitian diagonalization is appropriate. There is no B^2
    diamagnetic term or coupling to other torsion-rotation levels here.
    """
    moment, ids = magnetic_matrix_at_mf(j, species, vectors, fs, mf, gtr, gp)
    hamiltonian = np.diag(np.asarray(energies_khz)[ids]) - (
        MU_N_HZ_PER_MG / 1000 * field_mg) * moment
    energy, mixing = np.linalg.eigh(hamiltonian)
    return energy, mixing, ids
