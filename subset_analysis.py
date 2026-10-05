"""Export and summarize baseline melanoma PBMC samples from miraclib patients."""

import csv
import sqlite3
from contextlib import closing
from pathlib import Path


ROOT_FOLDER = Path(__file__).resolve().parent
DATABASE_FILE = ROOT_FOLDER / "cell_count.db"
OUTPUT_FILE = ROOT_FOLDER / "baseline_miraclib_melanoma_pbmc.csv"
SUMMARY_FILE = ROOT_FOLDER / "baseline_miraclib_melanoma_pbmc_summary.csv"
AVERAGE_FILE = ROOT_FOLDER / "melanoma_male_responder_baseline_b_cell_average.csv"

# All export and summary queries use the same baseline cohort definition.
ELIGIBLE_SAMPLES = """
SELECT
    sample,
    project,
    subject,
    response,
    sex,
    condition,
    treatment,
    sample_type,
    time_from_treatment_start AS timepoint,
    b_cell,
    cd8_t_cell,
    cd4_t_cell,
    nk_cell,
    monocyte
FROM cell_samples
WHERE condition = 'melanoma'
  AND sample_type = 'PBMC'
  AND treatment = 'miraclib'
  AND time_from_treatment_start = 0
"""

SUBSET_QUERY = f"""
WITH eligible_samples AS ({ELIGIBLE_SAMPLES}),
population_counts AS (
    SELECT sample, 'b_cell' AS population, b_cell AS cell_count
    FROM eligible_samples
    UNION ALL
    SELECT sample, 'cd8_t_cell', cd8_t_cell FROM eligible_samples
    UNION ALL
    SELECT sample, 'cd4_t_cell', cd4_t_cell FROM eligible_samples
    UNION ALL
    SELECT sample, 'nk_cell', nk_cell FROM eligible_samples
    UNION ALL
    SELECT sample, 'monocyte', monocyte FROM eligible_samples
)
SELECT
    samples.sample,
    samples.project,
    samples.subject,
    samples.response,
    samples.sex,
    samples.condition,
    samples.treatment,
    samples.sample_type,
    samples.timepoint,
    counts.population,
    counts.cell_count,
    totals.total_cell_count,
    100.0 * counts.cell_count / NULLIF(totals.total_cell_count, 0)
        AS percentage
FROM population_counts AS counts
JOIN eligible_samples AS samples ON samples.sample = counts.sample
JOIN (
    SELECT sample, SUM(cell_count) AS total_cell_count
    FROM population_counts
    GROUP BY sample
) AS totals ON totals.sample = counts.sample
ORDER BY samples.subject, samples.sample, counts.population
"""

# These are the three requested counts, calculated directly by SQLite.
SUMMARY_QUERIES = (
    (
        "samples_by_project",
        f"""
        WITH eligible_samples AS ({ELIGIBLE_SAMPLES})
        SELECT project AS group_name, COUNT(DISTINCT sample) AS count
        FROM eligible_samples
        GROUP BY project
        ORDER BY project
        """,
    ),
    (
        "subjects_by_response",
        f"""
        WITH eligible_samples AS ({ELIGIBLE_SAMPLES})
        SELECT response AS group_name, COUNT(DISTINCT subject) AS count
        FROM eligible_samples
        GROUP BY response
        ORDER BY response
        """,
    ),
    (
        "subjects_by_sex",
        f"""
        WITH eligible_samples AS ({ELIGIBLE_SAMPLES})
        SELECT sex AS group_name, COUNT(DISTINCT subject) AS count
        FROM eligible_samples
        GROUP BY sex
        ORDER BY sex
        """,
    ),
)

# This is intentionally a separate query from the PBMC/miraclib subset:
# include all melanoma sample and treatment types, limiting only sex, response,
# and baseline timepoint as requested.
MALE_RESPONDER_B_CELL_AVERAGE_QUERY = """
SELECT
    COUNT(*) AS sample_count,
    ROUND(AVG(b_cell), 2) AS average_b_cell_count
FROM cell_samples
WHERE condition = 'melanoma'
  AND sex = 'M'
  AND response = 'yes'
  AND time_from_treatment_start = 0
"""

OUTPUT_COLUMNS = (
    "sample",
    "project",
    "subject",
    "response",
    "sex",
    "condition",
    "treatment",
    "sample_type",
    "timepoint",
    "population",
    "cell_count",
    "total_cell_count",
    "percentage",
)


def export_baseline_subset():
    """Export sample-level data and SQL count results, verifying both outputs."""
    if not DATABASE_FILE.exists():
        raise FileNotFoundError(
            f"Database not found: {DATABASE_FILE}. Run load_data.py first."
        )

    # Write the five population rows for every eligible sample.
    with closing(sqlite3.connect(DATABASE_FILE)) as connection:
        cursor = connection.execute(SUBSET_QUERY)
        with OUTPUT_FILE.open("w", newline="", encoding="utf-8") as output_csv:
            writer = csv.writer(output_csv, lineterminator="\n")
            writer.writerow(OUTPUT_COLUMNS)
            writer.writerows(cursor)

        # Save each aggregate query's returned rows in one reviewable CSV.
        summary_rows = []
        for summary_name, query in SUMMARY_QUERIES:
            for group_name, count in connection.execute(query):
                summary_rows.append((summary_name, group_name, count))
        average_result = connection.execute(
            MALE_RESPONDER_B_CELL_AVERAGE_QUERY
        ).fetchone()

    if not summary_rows:
        raise ValueError(
            "No baseline melanoma PBMC samples from miraclib-treated patients "
            "were found."
        )

    with SUMMARY_FILE.open("w", newline="", encoding="utf-8") as summary_csv:
        writer = csv.writer(summary_csv, lineterminator="\n")
        writer.writerow(("summary", "group", "count"))
        writer.writerows(summary_rows)

    # Independently count distinct IDs in the detailed export and compare with
    # the SQL aggregate output, making it easy to verify both files agree.
    sample_groups = {}
    response_groups = {}
    sex_groups = {}
    with OUTPUT_FILE.open(newline="", encoding="utf-8") as subset_csv:
        for row in csv.DictReader(subset_csv):
            sample_groups.setdefault(row["project"], set()).add(row["sample"])
            response_groups.setdefault(row["response"], set()).add(row["subject"])
            sex_groups.setdefault(row["sex"], set()).add(row["subject"])

    exported_counts = {
        ("samples_by_project", group): len(ids)
        for group, ids in sample_groups.items()
    }
    exported_counts.update(
        {
            ("subjects_by_response", group): len(ids)
            for group, ids in response_groups.items()
        }
    )
    exported_counts.update(
        {
            ("subjects_by_sex", group): len(ids)
            for group, ids in sex_groups.items()
        }
    )
    sql_counts = {
        (summary, group): int(count)
        for summary, group, count in summary_rows
    }
    if sql_counts != exported_counts:
        raise ValueError(
            "SQL summary counts do not match distinct counts from the subset "
            f"CSV: SQL={sql_counts}, CSV={exported_counts}."
        )

    # Store the separate all-sample-types/treatments answer independently of
    # the baseline miraclib PBMC subset counts.
    average_sample_count, average_b_cell_count = average_result
    if average_sample_count == 0 or average_b_cell_count is None:
        raise ValueError(
            "No melanoma male responder samples were found at baseline."
        )
    with AVERAGE_FILE.open("w", newline="", encoding="utf-8") as average_csv:
        writer = csv.writer(average_csv, lineterminator="\n")
        writer.writerow(
            (
                "condition",
                "sex",
                "response",
                "timepoint",
                "sample_count",
                "average_b_cell_count",
            )
        )
        writer.writerow(
            (
                "melanoma",
                "M",
                "yes",
                0,
                average_sample_count,
                f"{average_b_cell_count:.2f}",
            )
        )

    print(f"Sample-level subset saved to: {OUTPUT_FILE}")
    print(f"SQL summary counts saved to: {SUMMARY_FILE}")
    print(f"Separate baseline B-cell average saved to: {AVERAGE_FILE}")
    print("SQL summary counts match the distinct counts in the subset CSV.")


def main():
    """Run the baseline subset export."""
    export_baseline_subset()


if __name__ == "__main__":
    main()
