# Immune Cell Response Analysis

This project imports `cell-count.csv` into SQLite, validates the source records,
calculates per-sample immune-cell relative frequencies, compares melanoma
responders with non-responders, and summarizes a baseline sample subset.

## Requirements

- Python 3.10 or later
- GNU Make
- `cell-count.csv` in the project root

## Setup

From the project root:

```sh
make setup
```

This installs dependencies from `requirements.txt`, including Streamlit for
the interactive dashboard.

## Run the full pipeline

```sh
make pipeline
```

The pipeline runs in order:

1. Inspect `cell-count.csv` and write `data_quality_report.csv`.
2. Stop if the source contains missing or malformed required values.
3. Load the validated data into `cell_count.db`.
4. Write Part 2 frequency results to `population_frequencies.csv`.
5. Run Part 3 statistical analysis and save its CSV, text interpretation, and
   boxplot.
6. Run Part 4 baseline subset queries and save the detailed sample export,
   SQL count summary, and separate male-responder B-cell average.

The pipeline does not impute or silently repair input values. Empty response
labels for healthy, untreated controls are treated as not applicable, matching
the loader's documented behavior.

## Start the dashboard

After setup and pipeline execution:

```sh
make dashboard
```

Streamlit prints a local URL (typically `http://localhost:8501`). The
interactive dashboard presents:

- **Part 2:** sample filters and population-frequency summaries from the
  pipeline-created SQLite database.
- **Part 3:** primary ANCOVA and supporting Welch results, their Holm-adjusted
  p-values, interpretation notes, and the population boxplot.
- **Part 4:** baseline samples by project, subjects by response and sex, and
  the average B-cell count for melanoma male responders at baseline.

The dashboard reads generated results and does not rebuild or modify the
database. Run `make pipeline` again after changing the input CSV.

## Statistical interpretation

Use Holm-adjusted p-values to judge statistical significance across each test
family; nominal p-values are exploratory when multiple comparisons are made.
The ANCOVA slope-interaction p-value is a diagnostic for the equal-slopes
assumption, not proof that all model assumptions hold. The Welch analyses
compare patient-level changes, while ANCOVA compares day-14 values adjusted
for baseline, so the two methods do not test identical contrasts.
