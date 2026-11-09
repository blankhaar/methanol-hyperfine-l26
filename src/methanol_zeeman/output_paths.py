"""Preflight generated products without modifying inputs or package data."""

from pathlib import Path


def _outside_package_data(path):
    resolved = path.resolve()
    package_parent = Path(__file__).resolve().parent.parent
    checkout=package_parent.parent
    if (checkout / "pyproject.toml").exists() and resolved.is_relative_to(checkout):
        raise ValueError("Choose an output outside the source checkout.")
    for package in ("methanol_hyperfine", "methanol_zeeman", "methanol_line_shift", "methanol_release"):
        protected = (package_parent / package).resolve()
        if resolved.is_relative_to(protected):
            raise ValueError("Choose an output path outside all immutable package data directories.")


def validate_output_file(path):
    """Reject an existing target or any destination under packaged input data.

    Resolve ancestor symlinks for the protected-data check. A dangling target
    symlink also counts as an existing target. Callers must still use exclusive
    creation (mode 'x'/'xb') to close the preflight-to-write race for files.
    """
    path = Path(path)
    _outside_package_data(path)
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Output already exists; choose a new file: {path}")
    for parent in path.parents:
        if (parent.exists() or parent.is_symlink()) and not parent.is_dir():
            raise NotADirectoryError(f"Output parent is not a directory: {parent}")
    return path


def validate_output_directory(path, filenames):
    """Validate every fixed product before computation; create nothing here."""
    path = Path(path)
    _outside_package_data(path)
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Choose a new output directory: {path}")
    for filename in filenames:
        validate_output_file(path / filename)
    return path
