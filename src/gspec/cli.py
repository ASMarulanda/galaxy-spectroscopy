"""Optional command-line interface.

The package's primary surface is the Python API (:func:`gspec.measure`). This
module exists so the pipeline can also be run as a program that writes one CSV —
useful for a batch job, for comparing against the validated baseline, and for
anyone who prefers files to objects.

The CSV schema is fixed and is *not* a contract the library owes anyone: it is
one serialisation of :meth:`gspec.results.SpectrumResult.row`. Values are written
with ``%.10g``, which preserves the physics; the input spectra are stored as
float32, so the 8th-10th printed digits are not data (DATA.md §5).

    gspec measure --cases data/cases --out out/line_measurements.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import variants
from .io import read_case
from .pipeline import measure

CSV_COLUMNS = ["case_id", "line", "redshift", "redshift_err", "flux", "flux_err", "flag"]

#: Kept for compatibility with the validated baseline's formatting.
FLOAT_FMT = "{:.10g}"


def _fmt(value):
    if isinstance(value, float):
        return FLOAT_FMT.format(value)
    return value


def measure_directory(cases_dir: str | Path, out_path: str | Path,
                      variant: str = variants.DEFAULT, verbose: bool = True) -> int:
    """Measure every ``*.fits`` in a directory and write the CSV. Returns the row count."""
    cases_dir = Path(cases_dir)
    files = sorted(cases_dir.glob("*.fits"))
    if not files:
        raise SystemExit(f"no case files under {cases_dir}")

    model = variants.get(variant)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_rows = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for path in files:
            case = read_case(path)
            result = measure(case, model=model)
            for line in result.lines:
                writer.writerow({k: _fmt(v) for k, v in result.row(line).items()})
                n_rows += 1
            if verbose:
                print(f"  {result.case_id}: z={result.redshift:.5f}  "
                      f"{result.model_selected}", flush=True)
    if verbose:
        print(f"wrote {out_path} ({n_rows} rows)")
    return n_rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gspec",
        description="Redshift and emission-line measurement for optical galaxy spectra.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    m = sub.add_parser("measure", help="measure a directory of case FITS files")
    m.add_argument("--cases", required=True, help="directory containing *.fits")
    m.add_argument("--out", required=True, help="CSV to write")
    m.add_argument("--variant", default=variants.DEFAULT, choices=sorted(variants.VARIANTS),
                   help="model configuration (default: m1a)")
    m.add_argument("--quiet", action="store_true")

    args = parser.parse_args(argv)
    if args.command == "measure":
        measure_directory(args.cases, args.out, args.variant, verbose=not args.quiet)
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
