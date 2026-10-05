"""
Inspect cell-count.csv for missing or malformed values without changing it.

Run `python data_quality.py` before `python load_data.py`. The script writes
one report row per issue to data_quality_report.csv; it never fills, removes,
or imputes any data.
"""

import csv
from collections import Counter
from pathlib import Path

from load_data import CSV_COLUMNS, CSV_FILE, INTEGER_COLUMNS


REPORT_FILE = Path(__file__).resolve().parent / "data_quality_report.csv"
SAMPLE_COLUMN_INDEX = CSV_COLUMNS.index("sample")


def add_issue(issues, record_number, sample, column, issue, value, details):
    """Add one consistent, easy-to-filter row to the report."""
    issues.append(
        {
            "csv_record": record_number,
            "sample": sample,
            "column": column,
            "issue": issue,
            "value": value,
            "details": details,
        }
    )


def inspect_csv():
    """Return a list of data issues and the number of CSV records inspected."""
    issues = []

    # csv.reader lets us compare the exact number of values in each record;
    # DictReader can otherwise hide whether a row had too many values.
    with CSV_FILE.open("r", newline="", encoding="utf-8") as csv_file:
        reader = csv.reader(csv_file)
        header = next(reader, None)

        if header is None:
            raise ValueError(f"The CSV file is empty: {CSV_FILE}")

        if header != list(CSV_COLUMNS):
            add_issue(
                issues,
                record_number=1,
                sample="",
                column="header",
                issue="unexpected_header",
                value=" | ".join(header),
                details=f"Expected: {' | '.join(CSV_COLUMNS)}",
            )

        # Header is record 1, so actual data records start at 2, matching the
        # line numbering used by load_data.py for this one-line-per-record CSV.
        records_inspected = 0
        for record_number, row in enumerate(reader, start=2):
            records_inspected += 1
            sample = (
                row[SAMPLE_COLUMN_INDEX].strip()
                if len(row) > SAMPLE_COLUMN_INDEX
                else ""
            )

            if len(row) < len(CSV_COLUMNS):
                # A short row has absent trailing fields. List each missing
                # expected column so it's clear what needs investigation.
                for column in CSV_COLUMNS[len(row) :]:
                    add_issue(
                        issues,
                        record_number,
                        sample,
                        column,
                        "missing_field",
                        "",
                        "This row has fewer values than the CSV header.",
                    )
            elif len(row) > len(CSV_COLUMNS):
                # Keep the extra raw values visible instead of silently
                # dropping them or guessing which column they belong to.
                for index, value in enumerate(
                    row[len(CSV_COLUMNS) :], start=len(CSV_COLUMNS) + 1
                ):
                    add_issue(
                        issues,
                        record_number,
                        sample,
                        f"extra_field_{index}",
                        "extra_field",
                        value,
                        "This value has no matching column in the CSV header.",
                    )

            # Check fields that are actually present. An empty string or
            # spaces-only string counts as missing, except for the documented
            # not-applicable response of untreated healthy control samples.
            for index, column in enumerate(CSV_COLUMNS[: len(row)]):
                value = row[index]
                if not value.strip():
                    response_not_applicable = (
                        column == "response"
                        and len(row) > max(
                            CSV_COLUMNS.index("condition"),
                            CSV_COLUMNS.index("treatment"),
                        )
                        and row[CSV_COLUMNS.index("condition")].strip() == "healthy"
                        and row[CSV_COLUMNS.index("treatment")].strip() == "none"
                    )
                    if response_not_applicable:
                        continue
                    add_issue(
                        issues,
                        record_number,
                        sample,
                        column,
                        "missing_value",
                        value,
                        "The CSV field is blank.",
                    )
                elif column in INTEGER_COLUMNS:
                    try:
                        int(value.strip())
                    except ValueError:
                        add_issue(
                            issues,
                            record_number,
                            sample,
                            column,
                            "invalid_integer",
                            value,
                            "Expected a whole-number value.",
                        )

    return issues, records_inspected


def write_report(issues):
    """Write issue details to a CSV file for sorting and filtering."""
    report_columns = (
        "csv_record",
        "sample",
        "column",
        "issue",
        "value",
        "details",
    )

    # newline="" lets the CSV library handle Windows line endings correctly.
    with REPORT_FILE.open("w", newline="", encoding="utf-8") as report_file:
        writer = csv.DictWriter(report_file, fieldnames=report_columns)
        writer.writeheader()
        writer.writerows(issues)


def main():
    """Inspect the original CSV and summarize issues without altering data."""
    issues, records_inspected = inspect_csv()
    write_report(issues)

    print(f"Inspected {records_inspected} CSV records.")
    if issues:
        issue_counts = Counter(issue["issue"] for issue in issues)
        print(f"Found {len(issues)} issues:")
        for issue_name, count in sorted(issue_counts.items()):
            print(f"  {issue_name}: {count}")
    else:
        print("No missing or malformed values found.")

    print(f"Detailed report written to: {REPORT_FILE}")
    print("No values were changed or imputed.")


if __name__ == "__main__":
    main()
