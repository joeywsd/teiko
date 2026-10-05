"""Run the complete data-quality, loading, analysis, and subset workflow."""

from contextlib import redirect_stdout
from pathlib import Path

import data_overview
import data_quality
import load_data
import statistical_analysis
import subset_analysis


ROOT_FOLDER = Path(__file__).resolve().parent
FREQUENCY_FILE = ROOT_FOLDER / "population_frequencies.csv"


def main():
    """Validate source rows, rebuild outputs, then run all analyses in order."""
    issues, records_inspected = data_quality.inspect_csv()
    data_quality.write_report(issues)
    print(
        f"Data quality: inspected {records_inspected} records; "
        f"found {len(issues)} issue(s)."
    )
    print(f"Data quality report: {data_quality.REPORT_FILE}")
    if issues:
        raise SystemExit(
            "Pipeline stopped: resolve the listed source CSV issues before "
            "loading the database."
        )

    # Loading only occurs after source checks pass; the loader replaces the
    # SQLite rows with the validated CSV contents.
    rows = load_data.read_csv_rows()
    load_data.create_database(rows)
    print(f"Loaded {len(rows)} samples into {load_data.DATABASE_FILE}.")

    # Persist Part 2's sample-level frequency query as CSV for reproducibility.
    with FREQUENCY_FILE.open("w", encoding="utf-8", newline="") as output_file:
        with redirect_stdout(output_file):
            data_overview.display_population_frequencies()
    print(f"Population frequencies saved to: {FREQUENCY_FILE}")

    # Run the existing Part 3 and Part 4 analysis scripts against the new DB.
    statistical_analysis.main()
    subset_analysis.main()


if __name__ == "__main__":
    main()
