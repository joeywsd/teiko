"""
Compare immune-cell population frequencies between treatment responders.

Run load_data.py first to create cell_count.db, then run this script.
The script keeps each patient's repeated measurements linked instead of
mistakenly treating three timepoints as three independent people.
"""

import csv
import io
import math
import sqlite3
from contextlib import closing, redirect_stdout
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf


ROOT_FOLDER = Path(__file__).resolve().parent
DATABASE_FILE = ROOT_FOLDER / "cell_count.db"
PLOT_FILE = ROOT_FOLDER / "responder_population_boxplots.png"
RESULTS_FILE = ROOT_FOLDER / "responder_statistics.csv"
SUMMARY_FILE = ROOT_FOLDER / "responder_analysis_summary.txt"

POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")
TIMEPOINTS = (0, 7, 14)
RESPONSE_LABELS = ("no", "yes")
SIGNIFICANCE_LEVEL = 0.05

# Keep percentages unrounded. Rounding is useful for display, but tests should
# use all the precision available from the original cell counts.
# The first CTE reshapes the five population columns into long-form rows,
# allowing each cell population to use the same downstream analysis.
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
    -- Compute the measured total for each sample before calculating shares.
    SELECT
        sample,
        population,
        population_count,
        SUM(population_count) OVER (PARTITION BY sample) AS total_count
    FROM population_counts
)
SELECT
    samples.subject,
    samples.response,
    samples.time_from_treatment_start AS timepoint,
    counts.sample,
    counts.population,
    counts.population_count AS cell_count,
    counts.total_count,
    -- Avoid dividing by zero; the loader reports the resulting missing value.
    100.0 * counts.population_count / NULLIF(counts.total_count, 0)
        AS percentage
FROM counts_with_totals AS counts
JOIN cell_samples AS samples ON samples.sample = counts.sample
-- Restrict the analysis to the treatment cohort being studied.
WHERE samples.condition = 'melanoma'
  AND samples.treatment = 'miraclib'
  AND samples.sample_type = 'PBMC'
ORDER BY samples.subject, counts.population, samples.time_from_treatment_start
"""


def load_frequency_data():
    """Load eligible patients' unrounded relative frequencies from SQLite."""
    # Fail early with an actionable message if data preparation has not run.
    if not DATABASE_FILE.exists():
        raise FileNotFoundError(
            f"Database not found: {DATABASE_FILE}. Run load_data.py first."
        )

    # closing() releases the connection after pandas finishes reading; the
    # sqlite3 transaction context alone would not close it.
    with closing(sqlite3.connect(DATABASE_FILE)) as connection:
        data = pd.read_sql_query(FREQUENCY_QUERY, connection)

    # Reject empty, invalid, or incomplete inputs before statistical analysis.
    if data.empty:
        raise ValueError(
            "No melanoma, miraclib, PBMC samples were found in the database."
        )
    if data["percentage"].isna().any():
        raise ValueError(
            "At least one eligible sample has a total cell count of zero; "
            "its relative frequencies cannot be calculated."
        )
    if not data["response"].isin(RESPONSE_LABELS).all():
        raise ValueError("Expected response values to be 'yes' or 'no'.")
    if set(data["response"].unique()) != set(RESPONSE_LABELS):
        raise ValueError(
            "Both responder and non-responder groups must be present."
        )
    if set(data["population"].unique()) != set(POPULATIONS):
        raise ValueError(
            f"Expected all five cell populations: {list(POPULATIONS)}."
        )

    return data


def validate_repeated_samples(data):
    """Check patient labels, uniqueness, and presence of all three visits."""
    # Response is a patient-level label. Conflicting labels would make it
    # unclear which response group that patient's repeated samples belong to.
    response_counts = data.groupby("subject")["response"].nunique()
    if (response_counts != 1).any():
        raise ValueError("A patient has conflicting response values.")

    # There should be exactly one value for each patient, population, and day.
    duplicate_rows = data.duplicated(
        subset=["subject", "population", "timepoint"]
    )
    if duplicate_rows.any():
        raise ValueError("Duplicate patient, population, and timepoint rows found.")

    # Require all three visits per cell population so changes can be calculated
    # directly, without imputing missing measurements.
    for population in POPULATIONS:
        # Inspect the visit schedule separately for each cell type.
        population_data = data[data["population"] == population]
        visits_per_patient = population_data.groupby("subject")[
            "timepoint"
        ].agg(lambda visits: set(visits))
        # Compare each patient's observed visits with the required schedule.
        incomplete_mask = visits_per_patient.map(
            lambda visits: visits != set(TIMEPOINTS)
        )
        incomplete = visits_per_patient[incomplete_mask]
        if not incomplete.empty:
            raise ValueError(
                f"{len(incomplete)} patients do not have all timepoints "
                f"for {population}."
            )


def plot_population_frequencies(data):
    """Save frequency boxplots grouped by response, population, and visit day."""
    # Allocate one panel per population; a sixth panel is left unused below.
    figure, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=True)
    axes = axes.flatten()

    # Color identifies visit day, while offsets cluster visits within response
    # groups without making their boxplots overlap.
    colors = {0: "#4C78A8", 7: "#F2A541", 14: "#59A14F"}
    offsets = {0: -0.22, 7: 0.0, 14: 0.22}

    for axis, population in zip(axes, POPULATIONS):
        # Restrict this panel to one cell population.
        population_data = data[data["population"] == population]
        box_values = []
        box_positions = []

        for response_index, response in enumerate(RESPONSE_LABELS):
            for timepoint in TIMEPOINTS:
                # Gather values for one response-group and visit combination.
                values = population_data.loc[
                    (population_data["response"] == response)
                    & (population_data["timepoint"] == timepoint),
                    "percentage",
                ].to_numpy()
                box_values.append(values)
                # Position the visit box beside its response-group partners.
                box_positions.append(
                    response_index * 1.5 + 1 + offsets[timepoint]
                )

        # Draw the six response-by-timepoint distributions for this panel.
        boxes = axis.boxplot(
            box_values,
            positions=box_positions,
            widths=0.19,
            patch_artist=True,
            showfliers=False,
        )
        # Keep timepoint colors consistent in every population panel.
        for box, timepoint in zip(
            boxes["boxes"], TIMEPOINTS * len(RESPONSE_LABELS)
        ):
            box.set_facecolor(colors[timepoint])
            box.set_alpha(0.75)

        axis.set_title(population)
        axis.set_xticks((1, 2.5), ("Non-responder", "Responder"))
        axis.set_ylabel("Cells in population (% of sample)")
        axis.grid(axis="y", alpha=0.25)

    # Hide the unused panel, then add a shared legend and figure title.
    axes[-1].axis("off")
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors[day], label=f"Day {day}")
        for day in TIMEPOINTS
    ]
    figure.legend(
        handles=legend_handles,
        title="Sample timepoint",
        loc="upper center",
        ncol=len(TIMEPOINTS),
    )
    figure.suptitle(
        "PBMC population frequencies by treatment response",
        fontsize=15,
        y=0.98,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    figure.savefig(PLOT_FILE, dpi=180)
    plt.close(figure)


def holm_adjust(p_values):
    """Apply Holm's step-down correction while preserving input order."""
    # Sort indices instead of values so adjusted results can be returned in
    # the same order as the tests that produced the original p-values.
    number_of_tests = len(p_values)
    sorted_indices = sorted(
        range(number_of_tests), key=lambda index: p_values[index]
    )
    adjusted = [0.0] * number_of_tests
    running_maximum = 0.0

    # Scale p-values using the number of hypotheses remaining at each rank.
    # The running maximum keeps adjusted p-values monotonic in sorted order.
    for rank, index in enumerate(sorted_indices):
        candidate = min(1.0, (number_of_tests - rank) * p_values[index])
        running_maximum = max(running_maximum, candidate)
        adjusted[index] = running_maximum

    return adjusted


def patient_timepoint_table(data, population):
    """Return one row per patient with response and three visit frequencies."""
    # Select one cell type, then spread its visit records across day columns.
    population_data = data[data["population"] == population]
    frequencies = population_data.pivot(
        index="subject", columns="timepoint", values="percentage"
    )
    # Response is patient-level, so retain one response label per patient.
    responses = population_data[["subject", "response"]].drop_duplicates()

    # Join by patient to preserve the relationship among repeated visits.
    patient_data = frequencies.join(responses.set_index("subject"))
    # Use readable day labels in the downstream model formulas.
    patient_data = patient_data.rename(
        columns={0: "day_0", 7: "day_7", 14: "day_14"}
    )
    # Encode response as 0/1 for use as a model predictor.
    patient_data["responder"] = (patient_data["response"] == "yes").astype(int)
    return patient_data


def run_ancova_tests(data):
    """Fit baseline-adjusted day-14 ANCOVA models for each population."""
    results = []

    for population in POPULATIONS:
        # Prepare repeated visit measurements for this cell type.
        patient_data = patient_timepoint_table(data, population)

        # Compare response groups at a common baseline frequency. HC3 gives
        # standard errors robust to unequal residual variances.
        model = smf.ols(
            "day_14 ~ responder + day_0", data=patient_data
        ).fit(cov_type="HC3")
        # Extract the adjusted responder difference and its 95% confidence
        # interval from the fitted model.
        estimate = model.params["responder"]
        interval = model.conf_int().loc["responder"]

        # Check the common-slope assumption by allowing the baseline-to-outcome
        # relationship to differ between response groups.
        slope_model = smf.ols(
            "day_14 ~ responder * day_0", data=patient_data
        ).fit(cov_type="HC3")

        # Keep the estimate, uncertainty, test result, and slope diagnostic
        # together for reporting and CSV export.
        results.append(
            {
                "analysis": "Day-14 ANCOVA (adjusted for day 0)",
                "population": population,
                "contrast": "Responder minus non-responder",
                "n_responder": int((patient_data["responder"] == 1).sum()),
                "n_non_responder": int((patient_data["responder"] == 0).sum()),
                "difference_percentage_points": estimate,
                "ci_lower": interval.iloc[0],
                "ci_upper": interval.iloc[1],
                "p_value": model.pvalues["responder"],
                "baseline_slope_interaction_p": slope_model.pvalues[
                    "responder:day_0"
                ],
            }
        )

    # Correct the five primary endpoint tests as one family.
    adjusted_values = holm_adjust([result["p_value"] for result in results])
    for result, adjusted_p in zip(results, adjusted_values):
        result["p_value_holm"] = adjusted_p

    return results


def welch_change_summary(responder_changes, non_responder_changes):
    """Compare mean changes using Welch's unequal-variance t-test."""
    # Summarize each response group's size, average change, and sample variance.
    responder_count = len(responder_changes)
    non_responder_count = len(non_responder_changes)
    responder_mean = responder_changes.mean()
    non_responder_mean = non_responder_changes.mean()
    responder_variance = responder_changes.var(ddof=1)
    non_responder_variance = non_responder_changes.var(ddof=1)
    # Explicit ddof=1 requests sample variances for the Welch standard error.
    standard_error_squared = (
        responder_variance / responder_count
        + non_responder_variance / non_responder_count
    )
    if standard_error_squared == 0:
        raise ValueError(
            "Cannot estimate a Welch test because both groups have zero "
            "variance in a change score."
        )

    # Welch-Satterthwaite degrees of freedom account for the two variances and
    # group sizes when building the confidence interval.
    degrees_of_freedom = standard_error_squared**2 / (
        (responder_variance / responder_count) ** 2 / (responder_count - 1)
        + (non_responder_variance / non_responder_count) ** 2
        / (non_responder_count - 1)
    )
    # Calculate the responder-minus-non-responder effect and its standard error.
    difference = responder_mean - non_responder_mean
    standard_error = standard_error_squared**0.5
    critical_value = stats.t.ppf(
        1 - SIGNIFICANCE_LEVEL / 2, degrees_of_freedom
    )
    test = stats.ttest_ind(
        responder_changes,
        non_responder_changes,
        equal_var=False,
        alternative="two-sided",
    )

    # Return statistics needed to report both the size and uncertainty of the
    # observed difference, plus the unadjusted test p-value.
    return {
        "n_responder": responder_count,
        "n_non_responder": non_responder_count,
        "mean_change_responder": responder_mean,
        "mean_change_non_responder": non_responder_mean,
        "difference_percentage_points": difference,
        "ci_lower": difference - critical_value * standard_error,
        "ci_upper": difference + critical_value * standard_error,
        "p_value": float(test.pvalue),
    }


def run_welch_change_tests(data):
    """Compare patient-level changes between responders and non-responders."""
    # Compare each adjacent-visit change and the overall day-0-to-day-14 change.
    intervals = (
        ("Day 0 to 7", "day_0", "day_7"),
        ("Day 7 to 14", "day_7", "day_14"),
        ("Cumulative day 0 to 14", "day_0", "day_14"),
    )
    results = []

    for population in POPULATIONS:
        # Keep all visit values for one patient together before taking
        # differences between visit days.
        patient_data = patient_timepoint_table(data, population)
        for interval_name, start_day, end_day in intervals:
            # Calculate each patient's change, then split changes by response.
            change = patient_data[end_day] - patient_data[start_day]
            result = welch_change_summary(
                change[patient_data["responder"] == 1],
                change[patient_data["responder"] == 0],
            )
            result.update(
                {
                    "analysis": "Welch test on patient-level change",
                    "population": population,
                    "contrast": interval_name,
                }
            )
            # Tag every result with its population and interval for the report.
            results.append(result)

    # Correct the 15 supporting tests as one family: five populations by
    # three time intervals.
    adjusted_values = holm_adjust([result["p_value"] for result in results])
    for result, adjusted_p in zip(results, adjusted_values):
        result["p_value_holm"] = adjusted_p

    return results


def display_results(ancova_results, change_results):
    """Report the slope diagnostic, primary ANCOVA, then supporting Welch tests."""
    columns = [
        "population",
        "contrast",
        "n_responder",
        "n_non_responder",
        "difference_percentage_points",
        "ci_lower",
        "ci_upper",
        "p_value",
        "p_value_holm",
    ]
    ancova_table = pd.DataFrame(ancova_results)
    change_table = pd.DataFrame(change_results)

    # Report the ANCOVA common-slope diagnostic before interpreting its
    # primary model. A non-significant interaction is not proof of equal slopes.
    print("\nANCOVA diagnostic: baseline-by-response interaction")
    for result in ancova_table.itertuples(index=False):
        interaction_p = result.baseline_slope_interaction_p
        if interaction_p < SIGNIFICANCE_LEVEL:
            interpretation = (
                "evidence of different baseline slopes; interpret the "
                "common-slope ANCOVA cautiously"
            )
        else:
            interpretation = (
                "no evidence of different baseline slopes detected by this test"
            )
        print(
            f"- {result.population}: interaction p={interaction_p:.4g}; "
            f"{interpretation}."
        )

    # The baseline-adjusted day-14 model is the primary analysis.
    print("\nPrimary analysis: day-14 ANCOVA adjusted for baseline")
    print(ancova_table[columns].to_string(index=False, float_format="%.4g"))

    # Use Holm-adjusted p-values to identify primary results meeting the
    # family-wise alpha threshold, while retaining effect sizes and intervals.
    significant = ancova_table[
        ancova_table["p_value_holm"] < SIGNIFICANCE_LEVEL
    ]
    print("\nPrimary findings after Holm correction (alpha = 0.05):")
    if significant.empty:
        print("No population met the corrected significance threshold.")
    else:
        for result in significant.itertuples(index=False):
            print(
                f"- {result.population}: adjusted day-14 difference "
                f"{result.difference_percentage_points:.3f} percentage points "
                f"(95% CI {result.ci_lower:.3f} to {result.ci_upper:.3f}; "
                f"Holm-adjusted p={result.p_value_holm:.4g})"
            )

    # Welch tests are supporting analyses and follow the primary ANCOVA output.
    print("\nSupporting analysis: Welch tests on patient-level changes")
    print(change_table[columns].to_string(index=False, float_format="%.4g"))

    # The five frequencies are compositional: they sum to 100% within sample,
    # so changes in one population are related to the others.
    print(
        "\nInterpretation note: population frequencies are related because the "
        "five percentages in each sample sum to 100%. These are associations "
        "and do not by themselves establish causation."
    )


def _print_saved_results_analysis(results_file):
    """Build the significance, diagnostics, and comparison report from CSV."""
    if not results_file.exists():
        raise FileNotFoundError(f"Statistical results file not found: {results_file}")

    # Read the saved results so this interpretation is based on the same
    # persistent output that can be shared and reviewed later.
    with results_file.open(newline="", encoding="utf-8") as csv_file:
        rows = list(csv.DictReader(csv_file))

    if not rows:
        raise ValueError(f"No statistical result rows found in {results_file}.")

    # Separate the primary endpoint tests from the supporting change tests.
    ancova_results = [
        row for row in rows
        if row["analysis"].startswith("Day-14 ANCOVA")
    ]
    welch_results = [
        row for row in rows
        if row["analysis"].startswith("Welch test")
    ]
    if not ancova_results or not welch_results:
        raise ValueError(
            "Results CSV must contain both Day-14 ANCOVA and Welch test rows."
        )

    def number(row, field):
        """Convert a required numeric CSV field and identify malformed values."""
        value = row.get(field, "")
        if value == "":
            raise ValueError(
                f"Missing {field!r} value for {row.get('population', 'unknown')} "
                f"({row.get('contrast', 'unknown contrast')})."
            )
        try:
            parsed = float(value)
        except ValueError as error:
            raise ValueError(
                f"Invalid {field!r} value {value!r} in {results_file}."
            ) from error
        if not math.isfinite(parsed):
            raise ValueError(
                f"Non-finite {field!r} value {value!r} in {results_file}."
            )
        return parsed

    alpha = SIGNIFICANCE_LEVEL
    ancova_raw = [row for row in ancova_results if number(row, "p_value") < alpha]
    ancova_adjusted = [
        row for row in ancova_results if number(row, "p_value_holm") < alpha
    ]
    welch_raw = [row for row in welch_results if number(row, "p_value") < alpha]
    welch_adjusted = [
        row for row in welch_results if number(row, "p_value_holm") < alpha
    ]

    def format_result(row):
        """Format a contrast with its effect, confidence interval, and p-values."""
        return (
            f"{row['population']} ({row['contrast']}): "
            f"difference {number(row, 'difference_percentage_points'):.3f} "
            f"percentage points (95% CI {number(row, 'ci_lower'):.3f} to "
            f"{number(row, 'ci_upper'):.3f}; raw p={number(row, 'p_value'):.4g}, "
            f"Holm-adjusted p={number(row, 'p_value_holm'):.4g}; "
            f"n={row['n_responder']} responders, "
            f"{row['n_non_responder']} non-responders)"
        )

    print(f"\nProgrammatic interpretation of {results_file.name} (alpha={alpha:g})")
    print(
        f"Primary ANCOVA: {len(ancova_raw)}/{len(ancova_results)} raw p-values "
        f"and {len(ancova_adjusted)}/{len(ancova_results)} Holm-adjusted "
        "p-values are below alpha."
    )
    if ancova_raw:
        print(
            "Nominal ANCOVA findings (not significant after correction "
            "unless stated):"
        )
        for row in ancova_raw:
            print(f"- {format_result(row)}")
    else:
        print("No ANCOVA result has a raw p-value below alpha.")

    if ancova_adjusted:
        print("ANCOVA findings significant after Holm correction:")
        for row in ancova_adjusted:
            print(f"- {format_result(row)}")
    else:
        print("No cell population has a Holm-adjusted significant ANCOVA result.")

    # The interaction p-value is the only ANCOVA assumption diagnostic stored
    # in this CSV; it cannot establish the full set of model assumptions.
    slope_flags = [
        (row, number(row, "baseline_slope_interaction_p"))
        for row in ancova_results
    ]
    slope_differences = [
        (row, p_value) for row, p_value in slope_flags if p_value < alpha
    ]
    print("\nAvailable ANCOVA equal-slopes diagnostic:")
    if slope_differences:
        for row, p_value in slope_differences:
            print(
                f"- {row['population']}: interaction p={p_value:.4g}; "
                "evidence of a baseline-slope difference, so interpret the "
                "common-slope ANCOVA cautiously."
            )
    else:
        print(
            "No interaction p-value is below alpha; this diagnostic detected "
            "no evidence of unequal baseline slopes, but does not prove the "
            "assumption holds."
        )
    print(
        "The CSV does not contain residual-level or patient-level data needed "
        "to assess other assumptions such as residual normality or influential "
        "observations."
    )

    print(
        f"\nSupporting Welch change tests: {len(welch_raw)}/{len(welch_results)} "
        f"raw p-values and {len(welch_adjusted)}/{len(welch_results)} "
        "Holm-adjusted p-values are below alpha."
    )
    if welch_raw:
        print(
            "Nominal Welch findings (not significant after correction "
            "unless stated):"
        )
        for row in welch_raw:
            print(f"- {format_result(row)}")
    else:
        print("No Welch result has a raw p-value below alpha.")
    if not welch_adjusted:
        print("No Welch change comparison is significant after Holm correction.")

    # Compare direction only for cell types with a nominal finding in both
    # analyses; emphasize these tests estimate different contrasts.
    ancova_nominal_by_population = {
        row["population"]: row for row in ancova_raw
    }
    welch_nominal_by_population = {}
    for row in welch_raw:
        welch_nominal_by_population.setdefault(row["population"], []).append(row)
    shared_populations = (
        ancova_nominal_by_population.keys() & welch_nominal_by_population.keys()
    )

    print("\nCross-method comparison:")
    if not shared_populations:
        print("No population has a nominal p<alpha finding in both methods.")
    else:
        for population in sorted(shared_populations):
            ancova_row = ancova_nominal_by_population[population]
            ancova_direction = number(
                ancova_row, "difference_percentage_points"
            )
            for welch_row in welch_nominal_by_population[population]:
                welch_direction = number(
                    welch_row, "difference_percentage_points"
                )
                direction_match = (
                    "same direction"
                    if ancova_direction * welch_direction > 0
                    else "opposite directions"
                    if ancova_direction * welch_direction < 0
                    else "one estimated difference is zero"
                )
                print(
                    f"- {population}: ANCOVA and {welch_row['contrast']} "
                    f"Welch result point in {direction_match}; neither "
                    "comparison is multiplicity-adjusted significant."
                )
    print(
        "ANCOVA estimates day-14 frequency adjusted for baseline; Welch tests "
        "compare changes over the stated intervals. They are related but not "
        "identical hypotheses. Treat unadjusted p<alpha findings as exploratory."
    )


def analyze_saved_results(results_file, summary_file):
    """Write the CSV interpretation report to a separate text file."""
    # Capture the report text rather than sending the interpretation to the
    # console, then persist it so it can be shared separately from the CSV.
    report = io.StringIO()
    with redirect_stdout(report):
        _print_saved_results_analysis(results_file)

    summary_file.write_text(report.getvalue(), encoding="utf-8")
    return summary_file


def main():
    """Check the cohort, create boxplots, and report statistical comparisons."""
    # Load the cohort and confirm its repeated measurements are analysis-ready.
    data = load_frequency_data()
    validate_repeated_samples(data)

    # Count each patient once for the cohort summary, rather than counting
    # their multiple population and visit rows.
    patient_counts = data[["subject", "response"]].drop_duplicates()
    response_counts = patient_counts["response"].value_counts()
    print(
        "Eligible patients: "
        f"{response_counts.get('no', 0)} non-responders, "
        f"{response_counts.get('yes', 0)} responders."
    )

    # Generate the plot, calculate both analysis families, and display results.
    plot_population_frequencies(data)
    ancova_results = run_ancova_tests(data)
    change_results = run_welch_change_tests(data)
    display_results(ancova_results, change_results)

    # Save full-precision results, including adjusted p-values and the ANCOVA
    # slope diagnostic, for later inspection.
    results = pd.DataFrame(ancova_results + change_results)
    results.to_csv(RESULTS_FILE, index=False)
    # Interpret the persistent CSV rather than relying only on in-memory tables.
    analyze_saved_results(RESULTS_FILE, SUMMARY_FILE)
    print(f"\nBoxplots saved to: {PLOT_FILE}")
    print(f"Statistical results saved to: {RESULTS_FILE}")
    print(f"Analysis summary saved to: {SUMMARY_FILE}")


if __name__ == "__main__":
    # Run the workflow only when this file is executed directly.
    main()