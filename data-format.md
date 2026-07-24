# Input data format

microplotr / microplotpy read three plain tab-separated (or comma-separated) files. No phyloseq objects, no QIIME2 artifacts, no BIOM — just flat tables you can open in Excel.

## 1. `feature_table.tsv` — abundance table

| Feature_ID | S1  | S2  | S3  | ... |
|------------|-----|-----|-----|-----|
| F1         | 327 | 57  | 12  | ... |
| F2         | 377 | 52  | 346 | ... |

- First column **must** be named `Feature_ID` (one row per ASV/OTU/species/gene).
- All remaining columns are samples; column name = `Sample_ID`.
- Values must be numeric and **non-negative** (raw counts or already-relative abundances both work; plots that need relative abundance compute it internally from raw counts).
- A feature or sample that is all zeros is allowed but triggers a warning (likely filtered out upstream, or a hint something didn't import right).

## 2. `taxonomy.tsv` — feature annotation

| Feature_ID | Domain   | Phylum         | Class               | Order            | Family            | Genus       |
|------------|----------|----------------|---------------------|------------------|-------------------|-------------|
| F1         | Bacteria | Proteobacteria | Gammaproteobacteria | Enterobacterales | Enterobacteriaceae| Escherichia |

- First column **must** be `Feature_ID` and its values must exactly match the `Feature_ID` column in `feature_table.tsv` (same set, no extras, no missing).
- Remaining columns are taxonomic ranks. Minimum required: `Phylum` and `Genus` (used by default plot settings — other ranks optional but recommended: `Domain`/`Kingdom`, `Class`, `Order`, `Family`, `Species`).
- Missing/unclassified values are allowed — leave blank or write `Unclassified`; plots label these explicitly rather than dropping them silently.

## 3. `metadata.tsv` — sample annotation

| Sample_ID | Group | Depth_m | pH  |
|-----------|-------|---------|-----|
| S1        | A     | 5       | 6.8 |

- First column **must** be `Sample_ID` and its values must exactly match the sample columns in `feature_table.tsv` (same set, no extras, no missing).
- Remaining columns are your variables of interest:
  - **Categorical grouping columns** (e.g. `Group`, `Treatment`, `Site`) — needed for barplot faceting, alpha/beta diversity group comparisons. Need ≥2 levels.
  - **Continuous gradient columns** (e.g. `Depth_m`, `pH`, `Time_days`) — needed for the ASV-level gradient plot. Must parse as numeric.

## Validator

Before any plot is drawn, both the R package (`mp_validate()`) and the Python package (`validate()`) run the same checks and report **errors** (block plotting) and **warnings** (plotting still allowed):

| Check | Level |
|---|---|
| Required columns present in each file | error |
| `Feature_ID` / `Sample_ID` unique within their file | error |
| Sample IDs match exactly between `feature_table.tsv` and `metadata.tsv` | error |
| Feature IDs match exactly between `feature_table.tsv` and `taxonomy.tsv` | error |
| Abundance values numeric and non-negative | error |
| All-zero feature or sample | warning |
| Missing taxonomy values | warning (recorded as "Unclassified") |
| Gradient column not numeric (when a gradient plot is requested) | error |
| Grouping column has <2 levels (when a group comparison plot is requested) | error |

Each finding reports which file, which column/row, and a plain-English message — e.g. *"Sample 'S9' found in metadata.tsv but missing from feature_table.tsv columns."*

## Example datasets

`inst/extdata/` (R) and `examples/` (Python) ship identical example sets:
- `example_valid/` — clean dataset, 8 samples / 20 features, used in all vignette/docs examples.
- `example_broken_id_mismatch/` — metadata references a sample not present in the feature table.
- `example_broken_negative_counts/` — one negative abundance value.
- `example_broken_missing_columns/` — `taxonomy.tsv` missing the `Phylum` column.
- `example_broken_duplicate_ids/` — duplicated `Feature_ID` row in `feature_table.tsv`.

Use these to see exactly what the validator reports before trying your own data.
