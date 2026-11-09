"""BELGI states and native eigenfunction calculation with explicit reference access.

The default compiles/runs bundled BELGI Fortran. Stored eigenfunctions are
retained only for an explicitly requested reference comparison.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import numpy as np

# Scientific inputs are ordinary package files, not a sibling private project.
PROJECT = Path(__file__).resolve().parent

@dataclass(frozen=True, order=True)
class State:
    species: str
    vt: int
    j: int
    k: int
    parity: str = ""

    @property
    def key(self):
        return f"{self.species}_{self.vt}_{self.j}_{self.k}_{self.parity or 'E'}"

    @classmethod
    def from_key(cls, key):
        species, vt, j, k, parity = key.split('_')
        state = cls(species, int(vt), int(j), int(k), '' if species == 'E' and parity == 'E' else parity)
        state.request()
        return state

    def request(self):
        if self.species not in ("A", "E") or abs(self.k) > self.j or not 0<=self.j<=46 or not 0<=self.vt<=2:
            raise ValueError(self)
        if self.species == "A" and (self.k < 0 or self.parity not in ("+", "-")):
            raise ValueError(self)
        if self.species=='A' and self.k==0 and self.parity!='+':raise ValueError(self)
        if self.species=='E' and self.parity:raise ValueError(self)
        return f"{int(self.species == 'E')} {self.vt} {self.j} {self.k} {int(self.parity == '+')}"

@dataclass
class Wavefunction:
    state: State
    energy: float
    coefficients: np.ndarray  # rows K=-J..J, columns n=-ktronc..ktronc

    @property
    def norm(self):
        return float(np.sum(self.coefficients**2))

    @property
    def ktronc(self):
        return (self.coefficients.shape[1]-1)//2

    def moments(self):
        j = self.state.j
        k = np.arange(-j, j+1)
        m = 3*np.arange(-self.ktronc, self.ktronc+1) + int(self.state.species == "E")
        p = self.coefficients
        jz2 = float(np.sum(p*p*k[:,None]**2))
        pajz = float(np.sum(p*p*k[:,None]*m[None,:]))
        ladder2 = np.sqrt((j-k[:-2])*(j+k[:-2]+1)*(j-k[:-2]-1)*(j+k[:-2]+2))
        jp2 = float(np.sum(p[2:]*p[:-2]*ladder2[:,None]))
        return {"pajz": pajz, "jz2": jz2, "jp2": jp2, "norm": self.norm}

def export_wavefunctions(states, ktronc=10, force=False, *, backend="native", cache_dir=None):
    """Compute with BELGI, or explicitly load the packaged regression reference."""
    if backend == 'native':
        from .native_belgi import native_wavefunctions
        return native_wavefunctions(states, ktronc, force, cache_dir)
    if backend != 'reference':
        raise ValueError("Choose backend='native' or 'reference'")
    if cache_dir is not None:
        raise ValueError('The reference provider does not use a generated cache')
    states = sorted(set(states))
    if not states or not 4 <= ktronc <= 20:
        raise ValueError("Nonempty states and 4<=ktronc<=20 required")
    for state in states:
        state.request()
    if force:
        raise ValueError("force=True requires the native BELGI provider")
    directory = PROJECT / "data" / "wavefunctions"
    metadata = directory / f"ktronc-{ktronc}.json"
    if not metadata.exists():
        raise ValueError(f"No frozen wavefunctions for ktronc={ktronc}")
    meta = json.loads(metadata.read_text())
    missing = [s.key for s in states if s.key not in meta["energies"]]
    if missing:
        raise ValueError("States unavailable in frozen wavefunctions: " + ", ".join(missing))
    archive = directory / f"ktronc-{ktronc}.npz"
    if hashlib.sha256(archive.read_bytes()).hexdigest() != meta["archive_sha256"]:
        raise RuntimeError("Frozen wavefunction checksum mismatch")
    with np.load(archive, allow_pickle=False) as data:
        return {s: Wavefunction(s, meta["energies"][s.key], data[s.key].copy()) for s in states}


def main(argv=None):
    """Write native eigenfunctions and their reproducible source/build identity."""
    import argparse
    from methanol_zeeman.output_paths import validate_output_directory
    from .native_belgi import configuration
    parser = argparse.ArgumentParser(description=__doc__)
    select = parser.add_mutually_exclusive_group()
    select.add_argument('--state', action='append', help='State key; repeat for multiple states')
    select.add_argument('--all-reference-states', action='store_true', help='Compute every state in the reference state list')
    parser.add_argument('--ktronc', type=int, default=10)
    parser.add_argument('--cache-dir', type=Path)
    parser.add_argument('--force', action='store_true', help='Run BELGI even when a matching generated cache exists')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    out = validate_output_directory(args.output, ('wavefunctions.npz', 'manifest.json'))
    if out.exists():
        parser.error('--output must be a fresh directory')
    if args.state:
        states = [State.from_key(s) for s in args.state]
    elif args.all_reference_states:
        metadata = PROJECT / 'data/wavefunctions' / f'ktronc-{args.ktronc}.json'
        states = [State.from_key(s) for s in json.loads(metadata.read_text())['energies']]
    else:
        from methanol_release.forward import load_problem
        states = list(load_problem(wavefunctions='reference').d.prepared)
    waves = export_wavefunctions(states, args.ktronc, args.force, cache_dir=args.cache_dir)
    key, config, _, _ = configuration()
    out.mkdir(parents=True)
    np.savez_compressed(out / 'wavefunctions.npz', **{s.key: w.coefficients for s, w in waves.items()})
    report = dict(config, configuration_sha256=key, requested_force_run=args.force,
                  archived_wavefunctions_used=False, ktronc=args.ktronc,
                  energies={s.key: w.energy for s, w in waves.items()},
                  norms={s.key: w.norm for s, w in waves.items()},
                  archive_sha256=hashlib.sha256((out/'wavefunctions.npz').read_bytes()).hexdigest())
    (out/'manifest.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(provider='native BELGI', states=len(waves), output=str(out),
                         configuration_sha256=key), indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
