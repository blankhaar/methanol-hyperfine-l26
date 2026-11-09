"""Pinned BELGI source/constant interface, adapted from the project resource.

Hamiltonian matrix elements remain in the genuine bundled BELGI Fortran code.
The only native-source build patch is the reviewed SETUP array-clear bound.
Hougen, Kleiner & Godefroid (1994), doi:10.1006/jmsp.1994.1047;
Xu et al. (2008), doi:10.1016/j.jms.2008.03.017.
"""
from pathlib import Path
import hashlib, math, re, subprocess, sys

XU_FINAL_CONSTANTS_SHA256 = "11353561e41ef266477e76898282d46b58f7fe92f228cd1f41061b9c95716f3d"

XU_FINAL_NONZERO_COUNT = 119

BELGI_VENDORED_SOURCE_TREE_SHA256 = (
    "6ed8a782fd0dd2cbbe7ef01ca285ae700f8a50e91a0ad2d944a8ca04321a0e8e"
)

CST_PARAMETER_ORDER = tuple(
    """
    OA B C DJ DJK DK ODELN ODELK HJ HJK
    HKJ HK OHJ OHJK OHK FV OFV GV ALV AK1
    AK2 AK3 AK4 AK5 AK6 AK7 DAB ODAB DABJ DABK
    OLV C1 C2 C3 DUMMY1 AK5J C2J GVJ AK2J C1J
    ANV BK2 AMV DAC DBC C4 ODELTA DELTA AK1J AK5K
    AK2K AK1K AK3J AK3K BK1 AK7J C11 ODAB6 AK6K AK7K
    AK66 AK6J AK3B AK4B C2K C1K C4J C4K DELTAB ODELTB
    C12 FPARA RHO V3 V6 V9 RHORHO LJ LK LJKKK
    LJJKK LJJJK OLKKK OLJKK OLJJK OLJJJ HABK HABJ HABJK GJ6
    GJ4K2 GJ2K4 GK6 FK6 AK6KK AK6JK AK6JJ AK7KK AK7JK AK7JJ
    AK3BK AMUA0 AMUA3 AMUB0 AMUB3 AMUC3 AGJ6 AGJ4K2 AGJ2K4 AGK6
    FJ6 FJ4K2 FJ2K4 AK7BC V9J V9K ODABK C6J4 C6J2K2 C6K4
    G4J4 G4J2K2 G4K4 GAJ4 GAJ2K2 GAK4 AK3BJ AK4BJ AK4BK ODABJ
    ODAC ODBC DELTAJ DELTAK DAGJ DAGK CABK CABJ CBCK CBCJ
    DG4J DG4K DACK DBCJ DACJ V12 V12J AK4C AK3C DAC6
    DBC6 GB AK8 AK9 AK8J AK8K AK9J AK9K DELTJK DELTKK
    ODAB9 C9BC C6BCJ C6BCK AG3BCJ AG3BCK G4BCJ G4BCK AK4BC AK3BC
    AK7AB BC4 CORREC AK76 OBC4 AK3BK6 AK4BK4 AK3CK4 AK4D AK3D
    AK4CK AK3CK AK3BK4 AK4BK6 AK4CK4 DACJJ DACJK DACKK DBCK DBCJK
    ODABKK C2JK C1JK C4JK C6BCJK G3BCJK G4BCJK C2KK C1KK C4KK
    C6BCKK G3BCKK G4BCKK DG3J DG3K
    """.split()
)

BELGI_IMPLEMENTATION_FILES = (
    "setup.f", "rottor.f", "htorse.f", "aset.f", "hrtset.f",
    "h.f", "encal.f", "vset.f", "vrtdi.f", "resum.f", "matout.f",
    "iposa.f", "ipose.f", "parity.f", "irep.f", "encal2.f", "vtors.f",
    "rule.f", "printener.f", "printderivl.f", "expectation.f", "factt.f",
    "ad3j.f", "inten.f", "pythag.f", "blockdata.f",
)

_SETUP_BOUNDS_BUG = "      max = 2*nn2+1"

_SETUP_BOUNDS_FIX = "      max = kdmx"

_ASSIGNMENT = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9]*)\s*=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?)\s*,?\s*$"
)

class HamiltonianError(RuntimeError):
    """Raised when the frozen-Hamiltonian build or validation fails."""

def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def parse_final_constants(path: str | Path, *, require_xu_final: bool = True) -> dict[str, str]:
    """Parse named Fortran assignments without relying on positional files.

    Returned values retain their decimal spelling.  The Fortran namelist read
    therefore performs the same conversion to double precision as the
    historical calculation.
    """
    source = Path(path)
    assignments: dict[str, str] = {}
    for number, line in enumerate(source.read_text(encoding="ascii").splitlines(), 1):
        if not line.strip():
            continue
        match = _ASSIGNMENT.fullmatch(line)
        if not match:
            raise HamiltonianError(f"invalid constants assignment at {source}:{number}: {line!r}")
        name, literal = match.groups()
        name = name.upper()
        if name not in CST_PARAMETER_ORDER:
            raise HamiltonianError(f"unknown BELGI constant {name!r} at {source}:{number}")
        if name in assignments:
            raise HamiltonianError(f"duplicate BELGI constant {name!r} at {source}:{number}")
        value = float(literal.replace("D", "E").replace("d", "e"))
        if not math.isfinite(value):
            raise HamiltonianError(f"non-finite BELGI constant {name!r}")
        assignments[name] = literal

    if require_xu_final:
        digest = sha256_file(source)
        if digest != XU_FINAL_CONSTANTS_SHA256:
            raise HamiltonianError(
                f"constants SHA-256 {digest} is not the audited Xu final vector "
                f"{XU_FINAL_CONSTANTS_SHA256}"
            )
        if len(assignments) != XU_FINAL_NONZERO_COUNT:
            raise HamiltonianError(
                f"Xu final vector needs {XU_FINAL_NONZERO_COUNT} named assignments; "
                f"found {len(assignments)}"
            )
    return assignments

def render_complete_cst_vector(assignments: dict[str, str]) -> str:
    """Serialize the canonical 205-slot vector for the prediction driver."""
    unknown = set(assignments) - set(CST_PARAMETER_ORDER)
    if unknown:
        raise HamiltonianError(f"unknown BELGI constants: {sorted(unknown)}")
    values = [assignments.get(name, "0.0D0") for name in CST_PARAMETER_ORDER]
    return "\n".join(("BELGI_FROZEN_CST_205", *values, ""))

def _belgi_sources(source_dir: Path) -> tuple[list[Path], list[Path]]:
    implementations = [source_dir / name for name in BELGI_IMPLEMENTATION_FILES]
    includes = sorted(source_dir.glob("common_*.f"))
    missing = [path for path in implementations if not path.is_file()]
    if missing:
        raise HamiltonianError(f"missing BELGI source files: {[path.name for path in missing]}")
    if not includes:
        raise HamiltonianError(f"no BELGI common include files found under {source_dir}")
    return implementations, includes

def belgi_source_inventory(source_dir: Path) -> dict[str, object]:
    """Checksum every vendored Fortran source and its deterministic inventory."""
    hashes = {
        path.name: sha256_file(path)
        for path in sorted(source_dir.glob("*.f"), key=lambda path: path.name)
    }
    payload = "".join(f"{digest}  {name}\n" for name, digest in hashes.items())
    tree_sha256 = hashlib.sha256(payload.encode("ascii")).hexdigest()
    return {
        "file_count": len(hashes),
        "file_sha256": hashes,
        "tree_sha256": tree_sha256,
        "audited_bundled_tree_sha256": BELGI_VENDORED_SOURCE_TREE_SHA256,
        "matches_audited_bundled_tree": tree_sha256 == BELGI_VENDORED_SOURCE_TREE_SHA256,
    }

def _compile_driver(
    project_root: Path,
    source_dir: Path,
    build_dir: Path,
    compiler: Path,
    *,
    driver_filename: str = "belgi_frozen_levels.f",
    executable_name: str = "belgi_frozen_levels",
) -> tuple[Path, list[str], str, dict[str, object]]:
    implementations, _ = _belgi_sources(source_dir)
    compiled_implementations: list[Path] = []
    source_patches: dict[str, object] = {}
    for source in implementations:
        if source.name != "setup.f":
            compiled_implementations.append(source)
            continue
        setup_text = source.read_text(encoding="ascii")
        if setup_text.count(_SETUP_BOUNDS_BUG) != 1:
            raise HamiltonianError(
                "setup.f does not contain the audited VECTPR bounds defect exactly once"
            )
        compiled_setup = build_dir / "setup.f"
        compiled_setup.write_text(
            setup_text.replace(_SETUP_BOUNDS_BUG, _SETUP_BOUNDS_FIX),
            encoding="ascii",
        )
        compiled_implementations.append(compiled_setup)
        source_patches[source.name] = {
            "reason": (
                "prevent SETUP from clearing 95 entries in VECTPR's "
                "93-element KDMX dimension and corrupting ENCAL2 state tracking"
            ),
            "original_sha256": sha256_file(source),
            "compiled_sha256": sha256_file(compiled_setup),
            "replacement": f"{_SETUP_BOUNDS_BUG.strip()} -> {_SETUP_BOUNDS_FIX.strip()}",
        }
    driver = project_root / "fortran" / driver_filename
    if not driver.is_file():
        raise HamiltonianError(f"prediction driver is missing: {driver}")
    executable = build_dir / executable_name
    if sys.platform == "darwin":
        lapack_link_flags = ["-Wl,-framework", "-Wl,Accelerate"]
    else:
        lapack_link_flags = ["-llapack", "-lblas"]
    command = [
        str(compiler), "-std=legacy", "-ffixed-line-length-none",
        "-fallow-argument-mismatch", "-O2", "-I", str(source_dir),
        str(driver), *(str(path) for path in compiled_implementations),
        *lapack_link_flags, "-o", str(executable),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise HamiltonianError(
            "BELGI driver compilation failed\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    version = subprocess.run(
        [str(compiler), "--version"], text=True, capture_output=True, check=True
    ).stdout.splitlines()[0]
    return executable, command, version, source_patches
