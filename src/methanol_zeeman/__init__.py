"""L26 methanol hyperfine and weak-field magnetic predictions.

Uses measured rotational magnetic terms and calculated remaining tensor terms.
"""

from .catalog import calculate_catalog, contract_magnetic, compute_transition

__all__ = ["calculate_catalog", "contract_magnetic", "compute_transition"]
