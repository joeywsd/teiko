"""
Display each sample's immune-cell counts as relative frequencies.

Run load_data.py first to create the SQLite database, then run this script.
"""

import csv
import sqlite3
import sys
from pathlib import Path


# Keep the database next to this script, even when the script is started from
# another working directory.
ROOT_FOLDER = Path(__file__).resolve().parent
DATABASE_FILE = ROOT_FOLDER / "cell_count.db"

# The source table has one sample per row, with one column per cell type.
# SQLite has no general UNPIVOT command, so UNION ALL turns those five columns
# into five population rows for each sample.
FREQUENCY_QUERY = """
WITH population_counts AS (
    SELECT sample, 'b_cell' AS population, b_cell AS population_count
    FROM cell_samples
    UNION ALL
    SELECT sample, 'cd8_t_cell', cd8_t_cell
    FROM cell_samples
    UNION ALL
    SELECT sample, 'cd4_t_cell', cd4_t_cell
    FROM cell_samples
    UNION ALL
    SELECT sample, 'nk_cell', nk_cell
    FROM cell_samples
    UNION ALL
    SELECT sample, 'monocyte', monocyte
    FROM cell_samples
),
counts_with_totals AS (
    SELECT
        sample,
        population,
        population_count,
        SUM(population_count) OVER (PARTITION BY sample) AS total_count
    FROM population_counts
)
SELECT
    sample,
    total_count,
    population,
    population_count AS count,
    ROUND(100.0 * population_count / NULLIF(total_count, 0), 2) AS percentage
FROM counts_with_totals
ORDER BY
    sample,
    CASE population
        WHEN 'b_cell' THEN 1
        WHEN 'cd8_t_cell' THEN 2
        WHEN 'cd4_t_cell' THEN 3
        WHEN 'nk_cell' THEN 4
        WHEN 'monocyte' THEN 5
    END
"""


def display_population_frequencies():
    """Print the Part 2 frequency summary as CSV-style output."""
    # Do not let SQLite silently create an empty database if loading has not
    # happened yet; tell the user how to create the database first.
    if not DATABASE_FILE.exists():
        raise FileNotFoundError(
            f"Database not found: {DATABASE_FILE}. Run load_data.py first."
        )

    connection = sqlite3.connect(DATABASE_FILE)

    try:
        cursor = connection.execute(FREQUENCY_QUERY)

        # There are five output rows per sample. CSV-style output is much
        # easier to scan or save than lining up tens of thousands of spaces.
        writer = csv.writer(sys.stdout, lineterminator="\n")
        writer.writerow([column[0] for column in cursor.description])
        writer.writerows(cursor)
    finally:
        # Close the database file even if the query or printing fails.
        connection.close()


def main():
    """Run the overview report when this file is executed directly."""
    display_population_frequencies()


if __name__ == "__main__":
    main()