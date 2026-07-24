# ASV/taxon abundance over a continuous gradient

`mp_asv_gradient_plot()` (R) -- `taxa_gradient_plot()` (Python) -- relative
abundance (%) of the top-N taxa plotted against a continuous metadata
gradient (depth, pH, time...). Standard way to show which taxa respond to
an environmental gradient (cf. `phyloseq::psmelt()` + `ggplot(aes(x =
gradient, y = Abundance, color = Taxon))`, common in depth-profile /
time-series microbiome papers).

## Required input

Standard tidy input. `metadata.tsv` must include the numeric gradient
column you plot (`gradient_var`) -- checked as numeric by the validator's
`gradient_column` option (see `data-format.md`).

## Key parameters

| Param | Default | Effect |
|---|---|---|
| `gradient_var` | *(required)* | Numeric metadata column for the x-axis (e.g. `"Depth_m"`, `"pH"`). |
| `rank` | `"Genus"` | Taxonomy level plotted. |
| `group_rank` | `"Phylum"` | Upper level for the same nested/bold-header legend convention as the barplot. `NULL`/`None` for a flat legend. |
| `top_n`, `min_rel_abund`, `min_prevalence`, `detection` | see barplot | Same taxa filtering/pooling logic as `mp_taxa_barplot()` -- taxa excluded either way land in `other_label`. |
| `smooth` | `TRUE`/`True` | Per-taxon trend line (R: `geom_smooth(method = "loess")`; Python: LOWESS via `statsmodels`, falls back to a low-degree polynomial fit if `statsmodels` isn't installed). |
| `facet` | `FALSE`/`False` | One panel per taxon (small multiples) instead of one overlaid panel with a color legend -- clearer when many taxa are shown at once. |

## Styling

No title -- axis titles, points, optional trend line, and (unless
faceted) a legend. Faceted mode drops the legend and uses bold panel
labels instead (consistent with the barplot/bubbleplot facet convention).

## Example

```python
from microplotpy import load
from microplotpy.plots import taxa_gradient_plot

data = load("feature_table.tsv", "taxonomy.tsv", "metadata.tsv")
fig1 = taxa_gradient_plot(data, gradient_var="Depth_m", rank="Genus", group_rank="Phylum", top_n=10)
fig1.savefig("gradient.png", dpi=300, bbox_inches="tight")

fig2 = taxa_gradient_plot(data, gradient_var="Depth_m", rank="Genus", top_n=6, facet=True)
fig2.savefig("gradient_facet.png", dpi=300, bbox_inches="tight")
```

## Literature

- `phyloseq::psmelt()` + ggplot gradient plots: standard in depth-profile / environmental-gradient amplicon studies.
- LOESS/LOWESS smoothing for abundance trends: Cleveland WS 1979, *JASA*, 74(368):829-836.
