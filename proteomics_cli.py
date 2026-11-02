"""Run the local table-analysis engine without starting a web server."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


def run_analysis(input_path, config_path, output_path):
    """Analyze into a new directory without overwriting any previous run."""
    source = Path(input_path).expanduser().resolve(strict=True)
    settings = Path(config_path).expanduser().resolve(strict=True)
    output = Path(output_path).expanduser()
    if not source.is_file() or source.suffix.lower() not in (".xlsx", ".xls", ".csv", ".tsv"):
        raise ValueError("Input must be an Excel, CSV or TSV file.")
    if not settings.is_file() or settings.stat().st_size > 1024 * 1024:
        raise ValueError("Analysis configuration must be a JSON file of at most 1 MB.")
    config = json.loads(settings.read_text("utf-8"))
    if not isinstance(config, dict):
        raise ValueError("Analysis configuration must be a JSON object.")
    if output.is_symlink() or output.exists():
        raise ValueError("Output directory already exists. Choose a new path for this run.")
    output.mkdir(parents=True, exist_ok=False, mode=0o700)
    # This entry reads local files and bundled plot assets.
    from proteomics import analyze_file
    result = analyze_file(source, config, output)
    summary = {"summary": result["summary"], "files": result["files"], "warnings": result["warnings"]}
    (output / "run_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", "utf-8")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="Devil's Proteomics table analysis")
    parser.add_argument("--input", required=True, help="Excel, CSV or TSV input table")
    parser.add_argument("--config", required=True, help="Analysis settings as a JSON object")
    parser.add_argument("--output", required=True, help="New output directory")
    args = parser.parse_args(argv)
    try:
        result = run_analysis(args.input, args.config, args.output)
    except (ValueError, OSError):
        print("Analysis failed. Check the input, configuration and output path. Existing output folders cannot be reused.", file=sys.stderr)
        return 1
    print(json.dumps({"status": "Complete", "total_rows": result["summary"]["total_rows"],
                      "significant_rows": result["summary"]["significant_rows"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
