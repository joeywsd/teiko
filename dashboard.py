"""Interactive Streamlit dashboard for analysis Parts 2 through 4."""

import sqlite3
from pathlib import Path

import pandas as pd
import streamlit as st


ROOT_FOLDER = Path(__file__).resolve().parent
DATABASE_FILE = ROOT_FOLDER / "cell_count.db"
FREQUENCY_FILE = ROOT_FOLDER / "population_frequencies.csv"
STATISTICS_FILE = ROOT_FOLDER / "responder_statistics.csv"
STATISTICS_SUMMARY_FILE = ROOT_FOLDER / "responder_analysis_summary.txt"
SUBSET_SUMMARY_FILE = ROOT_FOLDER / "baseline_miraclib_melanoma_pbmc_summary.csv"
B_CELL_AVERAGE_FILE = (
    ROOT_FOLDER / "melanoma_male_responder_baseline_b_cell_average.csv"
)
PLOT_FILE = ROOT_FOLDER / "responder_population_boxplots.png"

POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")


@st.cache_data
def get_filter_options():
    """Read available Part 2 filter values from the pipeline database."""
    with sqlite3.connect(DATABASE_FILE) as connection:
        values = {}
        for column in ("project", "condition", "sample_type", "treatment"):
            values[column] = [
                row[0]
                for row in connection.execute(
                    f"SELECT DISTINCT {column} FROM cell_samples "
                    f"WHERE {column} IS NOT NULL ORDER BY {column}"
                )
            ]
        values["timepoint"] = [
            row[0]
            for row in connection.execute(
                "SELECT DISTINCT time_from_treatment_start FROM cell_samples "
                "ORDER BY time_from_treatment_start"
            )
        ]
    return values


@st.cache_data
def load_population_frequencies(filters):
    """Query sample-level population counts and relative frequencies."""
    population_selects = " UNION ALL ".join(
        f"SELECT sample, '{population}' AS population, {population} AS cell_count "
        "FROM cell_samples"
        for population in POPULATIONS
    )
    conditions = []
    parameters = []
    filter_columns = {
        "project": "project",
        "condition": "condition",
        "sample_type": "sample_type",
        "treatment": "treatment",
        "timepoint": "time_from_treatment_start",
    }
    for key, column in filter_columns.items():
        selected = filters[key]
        if selected:
            placeholders = ", ".join("?" for _ in selected)
            conditions.append(f"{column} IN ({placeholders})")
            parameters.extend(selected)
    where_clause = " AND ".join(conditions) if conditions else "1 = 1"
    query = f"""
        WITH population_counts AS ({population_selects}),
        totals AS (
            SELECT sample, SUM(cell_count) AS total_count
            FROM population_counts
            GROUP BY sample
        )
        SELECT
            samples.sample,
            samples.subject,
            samples.project,
            samples.condition,
            samples.treatment,
            samples.sample_type,
            samples.time_from_treatment_start AS timepoint,
            counts.population,
            counts.cell_count,
            totals.total_count,
            100.0 * counts.cell_count / NULLIF(totals.total_count, 0)
                AS percentage
        FROM population_counts AS counts
        JOIN totals ON totals.sample = counts.sample
        JOIN cell_samples AS samples ON samples.sample = counts.sample
        WHERE {where_clause}
        ORDER BY samples.sample, counts.population
    """
    with sqlite3.connect(DATABASE_FILE) as connection:
        return pd.read_sql_query(query, connection, params=parameters)


def show_part_2():
    """Display interactive population-frequency summaries from SQLite."""
    st.header("Part 2 · Cell population frequencies")
    st.write(
        "Explore measured cell counts and relative frequencies. Filters query "
        "the pipeline-created SQLite database."
    )
    options = get_filter_options()
    selections = {}

    with st.sidebar:
        st.subheader("Part 2 filters")
        for column in ("project", "condition", "sample_type", "treatment"):
            selections[column] = st.multiselect(
                column.replace("_", " ").title(),
                options[column],
                default=options[column],
            )
        selections["timepoint"] = st.multiselect(
            "Timepoint (days)",
            options["timepoint"],
            default=options["timepoint"],
        )

    data = load_population_frequencies(selections)
    if data.empty:
        st.info("No samples match the selected filters.")
        return

    sample_count = data["sample"].nunique()
    subject_count = data["subject"].nunique()
    first_metric, second_metric = st.columns(2)
    first_metric.metric("Samples", f"{sample_count:,}")
    second_metric.metric("Subjects", f"{subject_count:,}")

    population_summary = (
        data.groupby("population", as_index=False)
        .agg(
            mean_percentage=("percentage", "mean"),
            mean_cell_count=("cell_count", "mean"),
        )
        .sort_values("population")
    )
    st.subheader("Mean population frequency across selected samples")
    st.bar_chart(
        population_summary.set_index("population")[["mean_percentage"]]
    )
    st.dataframe(
        population_summary.rename(
            columns={
                "population": "Population",
                "mean_percentage": "Mean frequency (%)",
                "mean_cell_count": "Mean cell count",
            }
        ),
        hide_index=True,
        width="stretch",
    )
    with st.expander("Sample-level population data"):
        st.dataframe(data, hide_index=True, width="stretch")

    if FREQUENCY_FILE.exists():
        st.caption(
            "The full Part 2 export is also available as "
            f"`{FREQUENCY_FILE.name}`."
        )


def show_part_3():
    """Display saved primary and supporting responder analyses."""
    st.header("Part 3 · Responder analysis")
    if not STATISTICS_FILE.exists():
        st.info("Run `make pipeline` to generate responder statistics.")
        return

    statistics = pd.read_csv(STATISTICS_FILE)
    ancova = statistics[
        statistics["analysis"].str.startswith("Day-14 ANCOVA")
    ].copy()
    welch = statistics[
        statistics["analysis"].str.startswith("Welch test")
    ].copy()
    if "p_value_holm" in ancova:
        ancova["holm_significant"] = ancova["p_value_holm"] < 0.05
    if "p_value_holm" in welch:
        welch["holm_significant"] = welch["p_value_holm"] < 0.05

    primary, supporting = st.tabs(["Primary · ANCOVA", "Supporting · Welch"])
    with primary:
        st.write(
            "Day-14 frequency difference adjusted for baseline. Holm-adjusted "
            "p-values account for the five primary population tests."
        )
        st.dataframe(ancova, hide_index=True, width="stretch")
    with supporting:
        st.write(
            "Welch tests compare patient-level change scores across the listed "
            "intervals. Holm correction covers all 15 comparisons."
        )
        st.dataframe(welch, hide_index=True, width="stretch")

    if STATISTICS_SUMMARY_FILE.exists():
        with st.expander("Interpretation and limitations", expanded=True):
            st.text(STATISTICS_SUMMARY_FILE.read_text(encoding="utf-8"))
    else:
        st.warning(
            "The summary interpretation file is missing. Run `make pipeline` "
            "to regenerate it."
        )
    if PLOT_FILE.exists():
        st.image(PLOT_FILE, caption="Population frequencies by response")


def show_part_4():
    """Display baseline subset group counts and the separate B-cell average."""
    st.header("Part 4 · Baseline subset analysis")
    if not SUBSET_SUMMARY_FILE.exists():
        st.info("Run `make pipeline` to generate baseline subset summaries.")
        return

    summary = pd.read_csv(SUBSET_SUMMARY_FILE)
    for summary_name, heading in (
        ("samples_by_project", "Baseline samples by project"),
        ("subjects_by_response", "Baseline subjects by response"),
        ("subjects_by_sex", "Baseline subjects by sex"),
    ):
        st.subheader(heading)
        group_data = summary.loc[
            summary["summary"] == summary_name, ["group", "count"]
        ].rename(columns={"group": "Group", "count": "Count"})
        st.dataframe(group_data, hide_index=True, width="stretch")

    if B_CELL_AVERAGE_FILE.exists():
        average = pd.read_csv(B_CELL_AVERAGE_FILE)
        if not average.empty:
            result = average.iloc[0]
            st.metric(
                "Mean baseline B cells · melanoma male responders",
                f"{result['average_b_cell_count']:,.2f}",
                help=(
                    "Includes all sample and treatment types; filters only "
                    "melanoma, male, responder, and timepoint 0."
                ),
            )
            st.caption(
                f"Calculated from {int(result['sample_count']):,} samples; "
                "B-cell counts, not percentages."
            )


def main():
    """Render the dashboard, with a helpful message before the pipeline runs."""
    st.set_page_config(
        page_title="Immune Cell Response Analysis",
        page_icon="🧬",
        layout="wide",
    )
    st.title("Immune Cell Response Analysis")
    st.caption("Parts 2–4 · melanoma immune-cell population analysis")

    if not DATABASE_FILE.exists():
        st.error(
            "The analysis database is missing. Run `make setup` and then "
            "`make pipeline` before opening the dashboard."
        )
        st.stop()

    show_part_2()
    st.divider()
    show_part_3()
    st.divider()
    show_part_4()


if __name__ == "__main__":
    main()
