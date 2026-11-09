"""L26 forward model: exact hyperfine diagonalization and experimental responses."""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.linalg import helmert
from scipy.special import voigt_profile
from methanol_hyperfine.belgi import State
from methanol_hyperfine.tensors import state_hamiltonian
from methanol_hyperfine.spin import diagonalize, transition_components, projection, symmetry_blocks
from methanol_hyperfine.lambdip import components, prediction
from methanol_hyperfine.widths import sigma_gamma
DATA = Path(__file__).parent / 'data/fit'

class Dataset:

    def hamiltonians(self, p):
        return {s: state_hamiltonian(prep, p) for s, prep in self.prepared.items()}

    def xu_features(self, h):
        eigen = {s: diagonalize(s.j, s.species, h[s])[0] for s in {g[q] for g in self.xu for q in ['u', 'l']}}
        out = []
        for g in self.xu:
            if self.xu_observable == 'assigned-intensity':
                c = transition_components(g['u'].j, g['l'].j, 'E', h[g['u']], h[g['l']], return_indices=True)
                features = []
                for ids in g['feature_components']:
                    pairs = {(u, l) for i in ids for u in g['assignments']['u'][i] for l in g['assignments']['l'][i]}
                    selected = np.array([(int(row[6]), int(row[7])) in pairs for row in c])
                    if not selected.any():
                        raise ValueError(f"No allowed transition for {g['id']} feature {ids}")
                    features.append(np.average(c[selected, 0], weights=c[selected, 1]))
                f = np.array(features)
                out.append(f - f.mean())
                continue
            shifts = []
            for i in range(4):
                shifts.append(np.mean(eigen[g['u']][g['assignments']['u'][i]]) - np.mean(eigen[g['l']][g['assignments']['l'][i]]))
            f = np.array([np.mean(np.array(shifts)[ids]) for ids in g['feature_components']])
            out.append(f - f.mean())
        return out

    def belov_features(self, h):
        coeff = {s: np.trace(b @ h[s]).real for s, b in self.belov_contractions.items()}
        predicted = [abs(2 * (g['u'].j + 0.5) * coeff[g['u']] + 2 * (g['l'].j + 0.5) * coeff[g['l']]) for g in self.belov]
        return (np.array(predicted), np.ones(len(predicted)))

    def profile_prediction(self, g, h, nuisance):
        c = transition_components(g['u'].j, g['l'].j, g['species'], h[g['u']], h[g['l']])
        c[:, 1] /= c[:, 1].sum()
        centre, logsigma, loggamma, doppler = nuisance
        sigma, gamma = (np.exp(logsigma), np.exp(loggamma))
        frequencies = c[:, 0]
        strength = c[:, 1]
        shifts = [0.0] if g['doppler_replicas'] == 1 else [-doppler / 2, doppler / 2]
        columns = [voigt_profile(g['x'][:, None] - centre - shift - frequencies[None, :], sigma, gamma) @ strength for shift in shifts]
        columns.extend([np.ones(len(g['x'])), g['x'] / max(np.ptp(g['x']), 1)])
        design = np.column_stack(columns)
        coeff = np.linalg.lstsq(design[g['mask']], g['y'][g['mask']], rcond=None)[0]
        coeff[:len(shifts)] = np.maximum(coeff[:len(shifts)], 0.0)
        return (design @ coeff, coeff, c)

class SelectiveLikelihood:
    """Exact existing Belov/Xu residuals and analytic eigenvalue derivatives.

    Fixed F/PI blocks are diagonalized in batches. Xu feature means and
    Helmert contrasts are unchanged; no effective constants become new data.
    """

    def __init__(self, dataset):
        if dataset.xu_observable != 'published-means':
            raise ValueError('Analytic route is for the published feature definition')
        self.d = dataset
        self.np = len(dataset.default)
        self.bmask = np.array([g['published_weight'] > 0 for g in dataset.belov])
        a = []
        b = []
        for g in dataset.belov:
            aa = np.zeros(self.np)
            bb = 0.0
            for q in ('u', 'l'):
                s = g[q]
                fixed, basis = dataset.prepared[s]
                projector = dataset.belov_contractions[s]
                factor = 2 * (s.j + 0.5)
                aa += factor * np.einsum('ij,pji->p', projector, basis).real
                bb += factor * np.trace(projector @ fixed).real
            a.append(aa)
            b.append(bb)
        self.ba = np.array(a)
        self.b0 = np.array(b)
        self.bobs = np.array([g['observed_splitting_khz'] for g in dataset.belov])
        packs = {}
        self.offsets = {}
        offset = 0
        for state in sorted({g[q] for g in dataset.xu for q in ('u', 'l')}):
            self.offsets[state] = offset
            fixed, basis = dataset.prepared[state]
            for u, _, _ in symmetry_blocks(state.j, state.species):
                size = u.shape[1]
                packs.setdefault(size, []).append((np.arange(offset, offset + size), u.conj().T @ fixed @ u, np.einsum('ai,pab,bj->pij', u.conj(), basis, u)))
                offset += size
        self.packs = {size: (np.array([v[0] for v in values]), np.array([v[1] for v in values]), np.array([v[2] for v in values])) for size, values in packs.items()}
        matrix = []
        obs = []
        self.xslices = []
        self.centred_slices = []
        centred = []
        start = 0
        cstart = 0
        for g in dataset.xu:
            sticks = np.zeros((4, offset))
            for i in range(4):
                for q, sign in [('u', 1), ('l', -1)]:
                    ids = np.array(g['assignments'][q][i]) + self.offsets[g[q]]
                    sticks[i, ids] += sign / len(ids)
            means = np.array([sticks[ids].mean(axis=0) for ids in g['feature_components']])
            means -= means.mean(axis=0)
            h = helmert(len(means))
            matrix.extend(h @ means)
            oo = np.array([v['observed_khz'] for v in g['observations']])
            obs.extend(h @ oo)
            centred.extend(means)
            self.xslices.append(slice(start, start + len(means) - 1))
            start += len(means) - 1
            self.centred_slices.append(slice(cstart, cstart + len(means)))
            cstart += len(means)
        self.xmat = np.array(matrix)
        self.xobs = np.array(obs)
        self.centred = np.array(centred)
        self.nenergy = offset

    def energies(self, p, derivative=False):
        e = np.zeros(self.nenergy)
        jac = np.zeros((self.nenergy, self.np)) if derivative else None
        for size, (indices, fixed, basis) in self.packs.items():
            h = fixed + np.einsum('p,bpij->bij', p, basis)
            if derivative:
                values, u = np.linalg.eigh(h)
                grad = np.einsum('bik,bpij,bjk->bkp', u.conj(), basis, u).real
                jac[indices.ravel()] = grad.reshape(-1, self.np)
            else:
                values = np.linalg.eigvalsh(h)
            e[indices.ravel()] = values.ravel()
        return (e, jac) if derivative else e

    def residual(self, p, family='both', derivative=False):
        residual = []
        jac = []
        if family in ('belov', 'both'):
            signed = self.b0 + self.ba @ p
            residual.append((np.abs(signed) - self.bobs)[self.bmask])
            if derivative:
                jac.append((np.sign(signed)[:, None] * self.ba)[self.bmask])
        if family in ('xu', 'both'):
            if derivative:
                e, de = self.energies(p, True)
                jac.append(self.xmat @ de)
            else:
                e = self.energies(p)
            residual.append(self.xmat @ e - self.xobs)
        if not residual:
            raise ValueError(family)
        r = np.concatenate(residual)
        return (r, np.vstack(jac)) if derivative else r

class TransferProblem:

    @staticmethod
    def decode(v):
        sigma, gamma = sigma_gamma(np.exp(v[1]), v[2])
        return np.array([v[0], np.log(sigma), np.log(gamma), v[3] if len(v) == 4 else 0.0])

    def one_profile(self, g, h, n):
        pred, amp, components = self.d.profile_prediction(g, h, n)
        weight = np.sqrt(min(1.0, 50 / np.count_nonzero(g['mask']))) / g['signal_sigma']
        return ((pred - g['y'])[g['mask']] * weight, pred, amp, components)

class DirectProblem:

    def unpack(self, z):
        p = self.anchor.copy()
        p[self.free] = z[:len(self.free)]
        return p

    def molecular_profiles(self, p, z):
        h = self.d.hamiltonians(p)
        sticks = [components(g, h) for g in self.traces]
        direct = [prediction(g, c, z[sl], not self.args.no_crossovers) for g, c, sl in zip(self.traces, sticks, self.specs)]
        glob = []
        if self.global_problem:
            zz = z[self.global_slice]
            gp = self.global_problem
            glob = [gp.one_profile(g, h, gp.decode(zz[spec['slice']])) for g, spec in zip(self.d.profiles, gp.specs)]
        return (sticks, direct, glob)

    def evaluate(self, z, jac=False):
        if self._cache is None or not np.array_equal(z, self._cache['z']):
            p = self.unpack(z)
            r, j = self.s.residual(p, 'both', True)
            sticks, direct, glob = self.molecular_profiles(p, z)
            blocks = [r[self.keep]] + [v[0] for v in direct + glob]
            lengths = np.cumsum([0] + [len(v) for v in blocks])
            self._cache = dict(z=z.copy(), p=p, sticks=sticks, direct=direct, glob=glob, residual=np.concatenate(blocks), lengths=lengths, table_jac=j[self.keep][:, self.free])
        c = self._cache
        if jac and 'jacobian' not in c:
            out = np.zeros((len(c['residual']), len(z)))
            nt = c['lengths'][1]
            out[:nt, :len(self.free)] = c['table_jac']
            for col, idx in enumerate(self.free):
                step = 3e-06 * max(1.0, abs(z[col]))
                values = []
                for sign in [1, -1]:
                    p = c['p'].copy()
                    p[idx] += sign * step
                    _, direct, glob = self.molecular_profiles(p, z)
                    values.append(np.concatenate([v[0] for v in direct + glob]) if direct or glob else np.array([]))
                out[nt:, col] = (values[0] - values[1]) / (2 * step)
            for i, (g, sl) in enumerate(zip(self.traces, self.specs)):
                rows = slice(c['lengths'][i + 1], c['lengths'][i + 2])
                v = z[sl]
                for k in range(3):
                    step = 3e-06 * max(1.0, abs(v[k]))
                    a = v.copy()
                    b = v.copy()
                    a[k] += step
                    b[k] -= step
                    out[rows, sl.start + k] = (prediction(g, c['sticks'][i], a, not self.args.no_crossovers)[0] - prediction(g, c['sticks'][i], b, not self.args.no_crossovers)[0]) / (2 * step)
            if self.global_problem:
                gp = self.global_problem
                zz = z[self.global_slice]
                h = self.d.hamiltonians(c['p'])
                for i, (g, spec) in enumerate(zip(self.d.profiles, gp.specs)):
                    sl = spec['slice']
                    v = zz[sl]
                    ridx = 1 + len(self.traces) + i
                    rows = slice(c['lengths'][ridx], c['lengths'][ridx + 1])
                    for k in range(len(v)):
                        step = 3e-06 * max(1.0, abs(v[k]))
                        a = v.copy()
                        b = v.copy()
                        a[k] += step
                        b[k] -= step
                        out[rows, self.global_slice.start + sl.start + k] = (gp.one_profile(g, h, gp.decode(a))[0] - gp.one_profile(g, h, gp.decode(b))[0]) / (2 * step)
            c['jacobian'] = out
        return c

    def residual(self, z):
        return self.evaluate(z)['residual']

    def jacobian(self, z):
        return self.evaluate(z, True)['jacobian']

def load_problem(data=DATA, *, wavefunctions="native"):
    """Load observations and build operators from native BELGI eigenfunctions.

    wavefunctions="reference" explicitly replays the archived operator arrays.
    """
    with np.load(data / 'dataset.npz', allow_pickle=False) as archive:
        arrays = {k: archive[k] for k in archive.files}

    def decode(v):
        if isinstance(v, list):
            return [decode(x) for x in v]
        if not isinstance(v, dict):
            return v
        if set(v) == {'array'}:
            return arrays[v['array']].copy()
        if set(v) == {'state'}:
            return State(*v['state'])
        if set(v) == {'slice'} and isinstance(v['slice'], list):
            return slice(*v['slice'])
        if set(v) == {'tuple'}:
            return tuple((decode(x) for x in v['tuple']))
        if set(v) == {'mapping'}:
            return {decode(k): decode(x) for k, x in v['mapping']}
        return {k: decode(x) for k, x in v.items()}
    values = decode(json.loads((data / 'dataset.json').read_text()))
    d = Dataset()
    d.__dict__.update(values['dataset'])
    if wavefunctions == 'native':
        from methanol_hyperfine.belgi import export_wavefunctions
        from .operators import prepare_masked_state
        waves = export_wavefunctions(d.prepared)
        rebuilt = {}
        for state, wave in waves.items():
            fixed, basis = prepare_masked_state(wave, d.frame)
            rebuilt[state] = (fixed, basis * d.coordinate_to_fraction[:, None, None])
        d.prepared = rebuilt
    elif wavefunctions != 'reference':
        raise ValueError('Choose native or reference wavefunctions')
    d.wavefunction_source = wavefunctions
    p = DirectProblem()
    p.__dict__.update(values['problem'])
    p.args = SimpleNamespace(**p.args)
    p.d = d
    p.s = SelectiveLikelihood(d)
    p._cache = None
    gp = TransferProblem()
    gp.d = d
    gp.specs = values['low_specs']
    p.global_problem = gp
    from .evaluation_cache import install_exact_cache
    install_exact_cache(p)
    return p

def reduced_transform(data=DATA):
    """Unscaled 138 x 127 tangent: four d, seven M, two DD, responses."""
    with np.load(data / 'basis.npz', allow_pickle=False) as b:
        transform = np.zeros((138, 127))
        transform[:4, :4] = np.eye(4)
        transform[4:22, 4:11] = b['M_basis']
        transform[22:, 11:] = np.eye(116)
    return transform
