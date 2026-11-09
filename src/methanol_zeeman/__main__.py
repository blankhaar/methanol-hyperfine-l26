"""Calculate magnetic components with native BELGI and the fixed L26 tensors."""

import argparse
import json
from .catalog import calculate_catalog


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Directory for generated CSV and provenance")
    parser.add_argument("--family", action="append", help="Family label, e.g. 6.7 or 12.2; repeatable")
    parser.add_argument("--rebuild-components", action="store_true", help="Compatibility option; components are always rebuilt")
    parser.add_argument("--reconstruct-hyperfine", action="store_true", help="Compatibility option; L26 tensors are always contracted with native waves")
    args = parser.parse_args(argv)
    result = calculate_catalog(args.output, args.family, args.rebuild_components, args.reconstruct_hyperfine)
    print(json.dumps({key: result[key] for key in ("components", "states", "families", "fit_performed", "eigenstates_recomputed", "hyperfine_recontracted")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
