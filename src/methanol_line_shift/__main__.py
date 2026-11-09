"""Compute a line shift from native BELGI/L26 or a supplied component CSV."""

import argparse
import json

from .estimate import estimate_shift, load_components
from .line_pull_surrogates import line_pull_Hz
from .profile_estimators import SHAPES
from methanol_zeeman.output_paths import validate_output_file


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default="12.2GHz", help="Alias or exact transition identifier")
    parser.add_argument("--components", help="CSV with offset_kHz,strength columns")
    parser.add_argument("--width", type=float, help="Component absorption FWHM in kHz")
    parser.add_argument("--response", choices=SHAPES, default="Gaussian")
    parser.add_argument("--surrogate", action="store_true", help="Use bounded frozen approximation (two targets; 30–2000 kHz)")
    parser.add_argument("--output", help="Optional JSON output path; result is also printed")
    parser.add_argument("--rebuild-catalog", metavar="DIRECTORY", help="Calculate with BELGI and generate all 136 response component patterns, then exit")
    args = parser.parse_args(argv)
    if args.rebuild_catalog:
        if args.width is not None or args.components or args.surrogate or args.output:
            parser.error("Use --rebuild-catalog alone; estimate shifts from its CSV in a separate invocation.")
        from .rebuild import rebuild_response_catalog
        print(json.dumps(rebuild_response_catalog(args.rebuild_catalog), indent=2))
        return 0
    if args.width is None:
        parser.error("--width is required for a line-shift estimate.")
    output = None
    if args.output:
        try:
            output = validate_output_file(args.output)
        except (ValueError, OSError) as error:
            parser.error(str(error))
    if args.surrogate:
        if args.components:
            parser.error("The frozen surrogate cannot be applied to a user component file.")
        result = dict(target=args.target, response=args.response, component_absorption_FWHM_kHz=args.width,
            shift_Hz=line_pull_Hz(args.width, args.target, args.response), method="frozen rational approximation")
    else:
        offsets, strengths = load_components(args.target, args.components)
        result = estimate_shift(offsets, strengths, args.width, args.response)
        result.update(target=args.target, method="full-real-axis symmetric-profile L2 fit",
            component_provider='user CSV' if args.components else 'native BELGI and L26 tensors')
    result["fit_scope"] = "Line-estimator nuisance parameters only; no hyperfine molecular refit"
    rendered = json.dumps(result, indent=2) + "\n"
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x") as handle:
            handle.write(rendered)
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
