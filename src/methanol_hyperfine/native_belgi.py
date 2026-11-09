"""Compile and run the genuine bundled BELGI prediction code.

Only locally generated, source/compiler-fingerprinted results are reused.
Packaged reference eigenfunctions are never a fallback for this provider.
"""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from .native_inputs import (
    BELGI_VENDORED_SOURCE_TREE_SHA256, HamiltonianError, _compile_driver,
    belgi_source_inventory, parse_final_constants, render_complete_cst_vector,
    sha256_file,
)

DATA = Path(__file__).resolve().parent / "data/belgi"
_MEMORY = {}


def configuration():
    """Validate scientific inputs before looking for any generated result."""
    inventory = belgi_source_inventory(DATA / "source")
    if (not inventory["matches_audited_bundled_tree"]
            or sha256_file(DATA / "source/SHA256SUMS") != BELGI_VENDORED_SOURCE_TREE_SHA256):
        raise HamiltonianError("Bundled BELGI source differs from its declared snapshot")
    constants = DATA / "inputs/xu2008-final-constants.txt"
    vector = render_complete_cst_vector(parse_final_constants(constants))
    requested = os.environ.get("METHANOL_FORTRAN_COMPILER", "gfortran")
    compiler = shutil.which(requested)
    if not compiler:
        raise HamiltonianError(
            "GNU Fortran is required to compute BELGI eigenfunctions. Install gfortran "
            "and LAPACK/BLAS (Apple Accelerate on macOS), or set "
            "METHANOL_FORTRAN_COMPILER. Archived waves are not used as a fallback."
        )
    compiler = Path(compiler).resolve()
    version = subprocess.run([str(compiler), "--version"], check=True, text=True,
                             capture_output=True, timeout=30).stdout.splitlines()[0]
    record = {
        "provider": "native BELGI", "source_tree_sha256": inventory["tree_sha256"],
        "source_files": inventory["file_count"], "constants_sha256": sha256_file(constants),
        "driver_sha256": sha256_file(DATA / "fortran/belgi_wavefunctions.f"),
        "adapter_sha256": sha256_file(Path(__file__)),
        "input_interface_sha256": sha256_file(Path(__file__).with_name("native_inputs.py")),
        "compiler": version, "compiler_binary_sha256": sha256_file(compiler),
        "platform": sys.platform, "vtmax": 2, "state_tracking_label": 1,
        "energy_unit": "cm^-1", "torsional_fourier_basis": "m = 3*n + sigma",
        "coefficient_convention": "native real coefficients; no renormalization",
        "native_limitation": "ENCAL exports zero for the ninth contracted torsional vector",
    }
    key = hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()
    return key, record, compiler, vector


def cache_root(explicit=None):
    root = Path(explicit or os.environ.get(
        "METHANOL_BELGI_CACHE", Path.home() / ".cache/methanol-hyperfine/belgi"
    )).expanduser().absolute()
    package = Path(__file__).resolve().parent.parent
    if (package.parent / 'pyproject.toml').is_file():
        package = package.parent
    if root.resolve().is_relative_to(package):
        raise ValueError("BELGI build/output cache must be outside the source package")
    if any(p.is_symlink() for p in (root, *root.parents)):
        raise ValueError("BELGI cache may not traverse a symlink")
    root.mkdir(parents=True, exist_ok=True)
    return root


@contextmanager
def _lock(path):
    with path.open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def _build(folder, compiler, config):
    build = folder / "build"
    receipt = build / "build.json"
    if receipt.exists():
        meta = json.loads(receipt.read_text())
        executable = build / "belgi_wavefunctions"
        if sha256_file(executable) != meta["executable_sha256"]:
            raise HamiltonianError("Generated BELGI executable checksum changed")
        return executable, meta
    with tempfile.TemporaryDirectory(prefix="build-", dir=folder) as work:
        work = Path(work)
        executable, command, version, patches = _compile_driver(
            DATA, DATA / "source", work, compiler,
            driver_filename="belgi_wavefunctions.f", executable_name="belgi_wavefunctions",
        )
        portable = ["{FORTRAN_COMPILER}"] + [
            str(x).replace(str(DATA), "{BELGI_DATA}").replace(str(work), "{BUILD}")
            for x in command[1:]
        ]
        meta = dict(config, compile_command=portable, source_patches=patches,
                    executable_sha256=sha256_file(executable))
        (work / "build.json").write_text(json.dumps(meta, indent=2) + "\n")
        # Rename the complete build while the process lock is held.
        os.rename(work, build)
    return build / "belgi_wavefunctions", meta


def _parse(stdout, states, ktronc):
    from .belgi import Wavefunction
    start, end = "BEGIN_WAVEFUNCTIONS\n", "END_WAVEFUNCTIONS"
    if stdout.count(start) != 1 or stdout.count(end) != 1:
        raise HamiltonianError("Missing/ambiguous native eigenfunction output markers")
    waves, seen = {}, {}
    current = None
    for line in stdout.split(start, 1)[1].split(end, 1)[0].splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "STATE":
            if len(fields) != 8:
                raise HamiltonianError("Malformed native STATE record")
            index = int(fields[1]) - 1
            if not 0 <= index < len(states):
                raise HamiltonianError("Native state index outside request")
            current = states[index]
            if current in waves or list(map(int, fields[2:7])) != list(map(int, current.request().split())):
                raise HamiltonianError("Native state identity does not match request")
            energy = float(fields[7])
            if not np.isfinite(energy):
                raise HamiltonianError("Nonfinite native energy")
            waves[current] = Wavefunction(current, energy, np.zeros((2*current.j+1, 2*ktronc+1)))
            seen[current] = set()
        else:
            if current is None or len(fields) != 3:
                raise HamiltonianError("Malformed native coefficient record")
            k, n = map(int, fields[:2]); value = float(fields[2])
            if not (-current.j <= k <= current.j and -ktronc <= n <= ktronc):
                raise HamiltonianError("Native coefficient index outside basis")
            if (k, n) in seen[current] or not np.isfinite(value):
                raise HamiltonianError("Duplicate/nonfinite native coefficient")
            seen[current].add((k, n))
            waves[current].coefficients[k+current.j, n+ktronc] = value
    if set(waves) != set(states):
        raise HamiltonianError("Native output omits a requested state")
    for state, wave in waves.items():
        if len(seen[state]) != wave.coefficients.size:
            raise HamiltonianError("Native output omits coefficients for " + state.key)
    return waves


def native_wavefunctions(states, ktronc=10, force=False, cache_dir=None):
    """Calculate requested eigenfunctions, reusing only validated generated caches."""
    from .belgi import Wavefunction
    states = sorted(set(states))
    if not states or not 4 <= ktronc <= 20:
        raise ValueError("Nonempty states and 4<=ktronc<=20 required")
    for state in states:
        state.request()
    key, config, compiler, vector = configuration()
    folder = cache_root(cache_dir) / key
    folder.mkdir(exist_ok=True)
    memory = _MEMORY.setdefault((str(folder), ktronc), {})
    if not force and all(s in memory for s in states):
        return {s: Wavefunction(s, memory[s].energy, memory[s].coefficients.copy()) for s in states}
    with _lock(folder / "native.lock"):
        executable, build = _build(folder, compiler, config)
        found = {} if force else {s: memory[s] for s in states if s in memory}
        if not force:
            for path in sorted(folder.glob(f"waves-k{ktronc}-*.json")):
                meta = json.loads(path.read_text())
                needed = [s for s in states if s not in found and s.key in meta["energies"]]
                if not needed:
                    continue
                archive = path.with_suffix(".npz")
                if meta["configuration_sha256"] != key or sha256_file(archive) != meta["archive_sha256"]:
                    raise HamiltonianError("Generated BELGI wavefunction cache checksum changed")
                with np.load(archive, allow_pickle=False) as data:
                    for s in needed:
                        c = data[s.key].copy()
                        if c.shape != (2*s.j+1, 2*ktronc+1) or not np.isfinite(c).all():
                            raise HamiltonianError("Invalid generated BELGI coefficient array")
                        found[s] = Wavefunction(s, meta["energies"][s.key], c)
        missing = [s for s in states if s not in found]
        if missing:
            request = "\n".join(s.request() for s in missing) + "\n"
            request_sha = hashlib.sha256(request.encode()).hexdigest()
            name = f"waves-k{ktronc}-{request_sha}"
            with tempfile.TemporaryDirectory(prefix="run-", dir=folder) as work:
                work = Path(work)
                req = work / "requests.txt"; req.write_text(request)
                command = [str(executable), str(max(s.j for s in missing)), str(ktronc), str(req)]
                run = subprocess.run(command, input=vector, text=True, capture_output=True,
                                     cwd=work, timeout=1800)
                if run.returncode:
                    raise HamiltonianError("BELGI execution failed: " + run.stderr[-4000:])
                fresh = _parse(run.stdout, missing, ktronc)
                archive = work / (name + ".npz")
                np.savez_compressed(archive, **{s.key: w.coefficients for s, w in fresh.items()})
                meta = dict(build, configuration_sha256=key, ktronc=ktronc,
                            request_sha256=request_sha, request_count=len(missing),
                            energies={s.key: w.energy for s, w in fresh.items()},
                            norms={s.key: w.norm for s, w in fresh.items()},
                            archive_sha256=sha256_file(archive),
                            command=["{BELGI_EXECUTABLE}", str(max(s.j for s in missing)), str(ktronc), "{REQUESTS}"],
                            archived_wavefunctions_used=False)
                record = work / (name + ".json")
                record.write_text(json.dumps(meta, indent=2) + "\n")
                os.replace(archive, folder / archive.name)
                os.replace(record, folder / record.name)
                (folder / (name + ".requests.txt")).write_text(request)
                (folder / (name + ".stderr.txt")).write_text(run.stderr)
                found.update(fresh)
        memory.update(found)
    return {s: Wavefunction(s, found[s].energy, found[s].coefficients.copy()) for s in states}
