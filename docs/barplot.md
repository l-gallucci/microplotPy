# Taxa barplot

`mp_taxa_barplot()` (R) / `taxa_barplot()` (Python) -- stacked relative-abundance
(%) barplot, one bar per sample, segments are the top-N taxa at a chosen rank.

## Required input

Standard tidy `feature_table.tsv` + `taxonomy.tsv` + `metadata.tsv` (see
`data-format.md`). `taxonomy.tsv` must include the rank you plot (default
`Genus`) and, if using the nested legend, the grouping rank above it
(default `Phylum`).

## Key parameters

| Param | Default | Effect |
|---|---|---|
| `rank` | `"Genus"` | Taxonomy level shown as bar segments. |
| `group_rank` | `"Phylum"` | Upper level used to order/color/group `rank`. Set `NULL`/`None` to disable. |
| `top_n` | 10 | Individual taxa shown (ranked by abundance among taxa passing the filters below); rest pooled into `other_label` ("Other", grey). `NULL`/`None` shows every taxon that passes the filters, uncapped. |
| `min_rel_abund` | 0 | Minimum mean relative abundance (%) across all samples for a taxon to be eligible; taxa below this go to `other_label` regardless of `top_n`. `0` disables. |
| `min_prevalence` | 0 | Minimum fraction of samples (`0-1`), or an absolute count if `> 1`, in which a taxon must be detected to be eligible. `0` disables. |
| `detection` | 0 | Relative-abundance (%) threshold above which a taxon counts as "detected" in a sample, for `min_prevalence`. |
| `nested_legend` | `TRUE`/`True` | Bold `group_rank` header above the first taxon of each group in the legend; requires `group_rank`. |
| `fix_taxonomy` | `TRUE`/`True` | Runs the tax-fix step first (see below) so no segment is left blank. |
| `facet_var` | `NULL`/`None` | Metadata column to facet by (e.g. treatment group). |
| `sample_order` | `NULL`/`None` | A metadata column name to sort samples by, or an explicit sample order. |

## Filtering

Which taxa count as "top" is decided in two independent stages: `min_rel_abund`
and `min_prevalence` first narrow the field to eligible taxa (a taxon must
pass *both* if both are set), then `top_n` picks the most abundant among
survivors. A taxon excluded at either stage -- failed a threshold, or simply
outside the top N of what passed -- ends up in `other_label` identically.
Set `top_n = NULL`/`None` to show every taxon that clears the thresholds,
uncapped. This is the standard core-microbiome filtering logic (cf.
`phyloseq::filter_taxa`, `microbiome::core`, Shade & Handelsman 2012,
*Environmental Microbiology* -- "beyond the Venn diagram").

## Missing/unknown taxonomy

Real 16S taxonomy tables routinely have blank, `NA`, or `uncultured` genus
calls. Left as-is, a genus-level barplot silently drops that abundance or
shows a meaningless blank segment. By default the barplot runs a "tax-fix"
step first: any unknown value at a rank is replaced by the **last known
ancestor rank**, so e.g. an unclassified genus inside family
*Lactobacillaceae* is labeled `Unclassified_Lactobacillaceae` rather than
`NA`. Cascading unknowns (genus *and* species both missing) anchor to the
same true ancestor rather than chaining onto each other's fabricated labels.
This mirrors the `tax_fix()` approach in the microViz R package (Barnett et
al. 2021, *Bioinformatics*).

## Nested/grouped legend

With `group_rank` set, taxa are ordered so every group's members stay
contiguous (both groups and within-group taxa ordered by descending
abundance), each group gets a base hue with its members shaded from that hue
(a family of related colors), and the legend shows the group name as a bold
header above the first member -- the convention used for microbiome barplots
in the R `ggnested` package. Set `nested_legend = FALSE` for a flat legend
(plain labels, one hue per taxon, no grouping).

## Styling

No title/subtitle/caption is ever drawn -- only axis titles, tick labels, and
the legend -- so the figure is publication-ready as returned.

## Example

```python
from microplotpy import load
from microplotpy.plots import taxa_barplot

data = load("feature_table.tsv", "taxonomy.tsv", "metadata.tsv")
fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=10)
fig.savefig("barplot.png", dpi=300, bbox_inches="tight")
```

## Literature

- Standard stacked taxa barplot: QIIME2 `taxa-bar-plot`, `phyloseq::plot_bar()`, MicrobiomeAnalyst.
- Nested/grouped legend: `ggnested` (gmteunisse/ggnested).
- Taxonomy fallback labeling: `tax_fix()` in `microViz` (Barnett DJM et al. 2021, *Bioinformatics*, 37(16)).
