from __future__ import annotations


class WignerBackendUnavailable(RuntimeError):
    pass


_initialized_max_two_j = 0
_initialized_wigner_type = 0


def _load_pywigxjpf():
    try:
        import pywigxjpf as wig
    except ModuleNotFoundError as exc:
        raise WignerBackendUnavailable(
            "pywigxjpf is not available in this interpreter. "
            "Activate python/.venv or install the dependency in your environment."
        ) from exc
    return wig


def initialize(max_two_j: int = 200, wigner_type: int = 9) -> None:
    global _initialized_max_two_j, _initialized_wigner_type
    wig = _load_pywigxjpf()
    if _initialized_max_two_j >= max_two_j and _initialized_wigner_type >= wigner_type:
        return
    wig.wig_table_init(max_two_j, wigner_type)
    wig.wig_temp_init(max_two_j)
    _initialized_max_two_j = max_two_j
    _initialized_wigner_type = wigner_type


def wig3jj(two_j1: int, two_j2: int, two_j3: int, two_m1: int, two_m2: int, two_m3: int) -> float:
    wig = _load_pywigxjpf()
    return float(wig.wig3jj(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3))


def wig6jj(two_j1: int, two_j2: int, two_j3: int, two_j4: int, two_j5: int, two_j6: int) -> float:
    wig = _load_pywigxjpf()
    return float(wig.wig6jj(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6))


def wig9jj(
    two_j1: int,
    two_j2: int,
    two_j3: int,
    two_j4: int,
    two_j5: int,
    two_j6: int,
    two_j7: int,
    two_j8: int,
    two_j9: int,
) -> float:
    wig = _load_pywigxjpf()
    return float(wig.wig9jj(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6, two_j7, two_j8, two_j9))

initialize()
