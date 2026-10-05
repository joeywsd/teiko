"""
Create a SQLite database from cell-count.csv.

The table has one row per sample. That shape matches questions like
"what are the cell counts in this sample?" without needing a table join.
"""

import csv
import sqlite3
from pathlib import Path


# Put both files next to this script. This works even if Python is launched
# from a different folder, because __file__ means "the location of this file."
ROOT_FOLDER = Path(__file__).resolve().parent
CSV_FILE = ROOT_FOLDER / "cell-count.csv"
DATABASE_FILE = ROOT_FOLDER / "cell_count.db"

# These names and their order mirror the CSV header row. Keeping the names in
# one place helps the CSV reader and SQL insert agree about the column order.
CSV_COLUMNS = (
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
)

# SQLite stores whole-number measurements as INTEGERs. The other fields are
# labels, so they are stored as TEXT.
INTEGER_COLUMNS = (
    "age",
    "time_from_treatment_start",
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
)


def read_csv_rows():
    """Read the CSV and turn each sample into values ready for SQLite."""
    rows = []

    # newline="" lets Python's CSV library handle line endings correctly.
    # utf-8 is a widely used text encoding and handles normal CSV text safely.
    with CSV_FILE.open("r", newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)

        # A header mismatch could put values under the wrong names, so stop
        # with a clear error instead of quietly loading incorrect data.
        if reader.fieldnames != list(CSV_COLUMNS):
            raise ValueError(
                "Unexpected CSV columns. "
                f"Expected {list(CSV_COLUMNS)}, found {reader.fieldnames}."
            )

        # DictReader gives us one dictionary per line, using the header names
        # as keys. The enumerate count starts at 2 because line 1 is the header.
        for line_number, csv_row in enumerate(reader, start=2):
            # DictReader uses a None key when a row has extra values, and a
            # None value when a row is missing a trailing field. Neither is a
            # legitimate blank response, so report malformed rows explicitly.
            if None in csv_row or any(
                csv_row.get(column) is None for column in CSV_COLUMNS
            ):
                raise ValueError(
                    f"CSV line {line_number} has a missing or extra value."
                )

            # Strip accidental spaces, then convert counts and ages from text
            # (how CSV represents everything) into whole numbers.
            clean_row = {
                column: csv_row[column].strip() for column in CSV_COLUMNS
            }

            # In this dataset, response is intentionally blank for healthy
            # patients who received no treatment. Preserve that meaning as
            # Python None, which SQLite stores as SQL NULL; don't guess yes/no.
            response_is_not_applicable = (
                clean_row["response"] == ""
                and clean_row["condition"] == "healthy"
                and clean_row["treatment"] == "none"
            )

            # Other blank fields are not explained by that rule and could
            # change the meaning of an analysis, so fail rather than inventing
            # a replacement such as zero or an average.
            blank_columns = [
                column
                for column in CSV_COLUMNS
                if clean_row[column] == ""
                and not (
                    column == "response" and response_is_not_applicable
                )
            ]
            if blank_columns:
                raise ValueError(
                    f"CSV line {line_number} has blank required value(s): "
                    f"{', '.join(blank_columns)}."
                )

            if response_is_not_applicable:
                clean_row["response"] = None

            for column in INTEGER_COLUMNS:
                try:
                    clean_row[column] = int(clean_row[column])
                except ValueError as error:
                    raise ValueError(
                        f"CSV line {line_number}: {column} must be a whole number."
                    ) from error

            # A tuple is like a fixed row of values, in the same order as the
            # SQL column list. Save it so we can insert all rows together.
            rows.append(tuple(clean_row[column] for column in CSV_COLUMNS))

    return rows


def create_database(rows):
    """Create the samples table and fill it with the CSV rows."""
    connection = sqlite3.connect(DATABASE_FILE)

    try:
        # The connection context makes the database changes one transaction:
        # either all rows are saved, or none are if something goes wrong.
        with connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS cell_samples (
                    project TEXT NOT NULL,
                    subject TEXT NOT NULL,
                    condition TEXT NOT NULL,
                    age INTEGER NOT NULL,
                    sex TEXT NOT NULL,
                    treatment TEXT NOT NULL,
                    response TEXT,
                    sample TEXT PRIMARY KEY,
                    sample_type TEXT NOT NULL,
                    time_from_treatment_start INTEGER NOT NULL,
                    b_cell INTEGER NOT NULL,
                    cd8_t_cell INTEGER NOT NULL,
                    cd4_t_cell INTEGER NOT NULL,
                    nk_cell INTEGER NOT NULL,
                    monocyte INTEGER NOT NULL
                )
                """
            )

            # Make the script safe to run again: replace the old imported rows
            # rather than failing because their sample IDs already exist.
            connection.execute("DELETE FROM cell_samples")

            # The ? marks are placeholders. SQLite safely puts each value into
            # its place, instead of treating data as part of the SQL command.
            column_names = ", ".join(CSV_COLUMNS)
            value_placeholders = ", ".join("?" for _ in CSV_COLUMNS)
            insert_sql = (
                f"INSERT INTO cell_samples ({column_names}) "
                f"VALUES ({value_placeholders})"
            )
            connection.executemany(insert_sql, rows)
    finally:
        # A with-connection transaction commits or rolls back, but does not
        # close the connection. Always close it when we are finished.
        connection.close()


def main():
    """Load the source rows into the SQLite database."""
    rows = read_csv_rows()
    create_database(rows)
    print(f"Loaded {len(rows)} samples into {DATABASE_FILE}")


# This makes `python load_data.py` run the loader, but importing this file
# from another Python script will not unexpectedly create or change a database.
if __name__ == "__main__":
    main()
