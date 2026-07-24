# Functional profile barplot/heatmap

`mp_function_barplot()` / `mp_function_heatmap()` (R) --
`function_barplot()` / `function_heatmap()` (Python) -- thin wrappers around
the taxa barplot/heatmap engine (`mp_taxa_barplot()`/`mp_taxa_heatmap()`,
see `barplot.md`/`heatmap.md`) for functional (not taxonomic) data: gene
family / KO / pathway abundance instead of taxon abundance. Same nested
legend, same tax-fix fallback for missing annotation, same filters -- just
different input names and defaults.

## Why a wrapper, not a new plot

A per-sample gene/KO count table (`Feature_ID` + per-sample counts) joined
to a functional annotation table (`Feature_ID` + hierarchy columns) is
shape-identical to `feature_table.tsv` + `taxonomy.tsv`. Duplicating the
whole barplot/heatmap engine for this would just be the same code twice.

## Required input

Two tables, not the 3-file taxonomy schema:

**`gene_count_table.tsv`** -- one row per gene/KO/feature, same shape as
`feature_table.tsv`: `Feature_ID` + one column per sample. Comes from any
read quantifier run against predicted genes/ORFs (e.g. featureCounts,
Salmon, CoverM) -- this library doesn't do quantification, only plotting.

**`function_annotation.tsv`** -- `Feature_ID` + functional hierarchy
columns, from whichever annotation tool you used:

- **eggNOG-mapper** (`*.emapper.annotations`): use its own column names
  directly -- `COG_category` (broad functional category, single letter,
  e.g. `"C"` = energy production; `"S"` = function unknown -- a real,
  distinct COG category, not the same thing as *no* hit), `KEGG_ko`,
  `Description`. Defaults: `rank = "KEGG_ko"`, `group_rank = "COG_category"`.
- **KofamScan** (`mapper.tsv`/detail-tsv): only gives a flat KO assignment,
  no broader category (`Feature_ID`/`KO`/`KO_definition`). Use
  `rank = "KO"`, `group_rank = NULL`/`None`, `nested_legend = FALSE`/`False`
  (nothing to nest under).

`metadata` is optional (only needed for `facet_var`/`sample_order` by a
metadata column) -- auto-stubbed from the sample columns if not supplied.

## Missing annotation

A gene with no COG hit at all (common -- eggNOG-mapper leaves `COG_category`
blank, not `"S"`, when nothing matched) is resolved the same way an
unclassified genus is on the taxonomic side: [mp_tax_fix()]/`tax_fix()`
labels it `"Unclassified"`. Internally this is wired via `tax_fix_ranks`
(new parameter added to the taxa engine alongside this feature) so the
fallback machinery works on arbitrary column names, not just the fixed 16S
rank set.

## Example

```python
import pandas as pd
from microplotpy.plots import function_barplot, function_heatmap

gene_counts = pd.read_csv("gene_count_table.tsv", sep="\t")
annotation = pd.read_csv("function_annotation.tsv", sep="\t")

# eggNOG-mapper shape
fig1 = function_barplot(gene_counts, annotation, rank="KEGG_ko", group_rank="COG_category")
fig2 = function_heatmap(gene_counts, annotation, rank="KEGG_ko")

# KofamScan shape (flat KO annotation, no broader category)
fig3 = function_barplot(gene_counts, annotation_kofam, rank="KO", group_rank=None, nested_legend=False)
```

## Literature

- eggNOG-mapper: Cantalapiedra CP et al. 2021, *Molecular Biology and Evolution*, 38(12):5825-5829.
- KofamScan/KofamKOALA: Aramaki T et al. 2020, *Bioinformatics*, 36(7):2251-2252.
- HUMAnN3 (an alternative, pre-computed pathway-abundance route rather than gene-level annotation): Beghini F et al. 2021, *eLife*, 10:e65088.
