# Taxa bubble (dot) plot

`mp_taxa_bubbleplot()` (R) / `taxa_bubbleplot()` (Python) -- one row per
taxon, one column per sample, dot size = relative abundance (%), dot color =
an upper taxonomy level.

## Required input

Same tidy input as the barplot.

## Key parameters

| Param | Default | Effect |
|---|---|---|
| `rank` | `"Genus"` | Taxonomy level shown as rows. |
| `group_rank` | `"Phylum"` | Upper level used to color and order `rank`. Set `NULL`/`None` for a single-color plot. |
| `top_n` | 25 | Taxa shown (ranked by total relative abundance among taxa passing the filters below). `NULL`/`None` shows every taxon that passes the filters, uncapped. |
| `min_rel_abund` | 0 | Minimum mean relative abundance (%) across all samples for a taxon to be shown. `0` disables. |
| `min_prevalence` | 0 | Minimum fraction of samples (`0-1`), or an absolute count if `> 1`, in which a taxon must be detected to be shown. `0` disables. |
| `detection` | 0 | Relative-abundance (%) threshold above which a taxon counts as "detected" in a sample, for `min_prevalence`. |
| `fix_taxonomy` | `TRUE`/`True` | Same tax-fix step as the barplot -- see `barplot.md`. |
| `facet_var` | `NULL`/`None` | Metadata column to facet by. |

## Filtering

Same two-stage logic as the barplot (`min_rel_abund`/`min_prevalence` narrow
candidates, `top_n` caps the survivors) -- see `barplot.md` for the full
rationale. Taxa that don't make the cut are simply not shown (no "Other"
row, since a bubbleplot isn't a composition that needs to sum to 100%).

## Complete grid

Every sample x taxon combination is drawn, including absent ones (dot size
0), rather than only drawing rows/columns where a taxon was detected --
so presence/absence and prevalence patterns read directly off the grid,
matching the `phyloseq::psmelt()` + ggplot bubble / MicrobiomeAnalyst
dot-plot convention. Rows are ordered grouped by `group_rank` block (blocks
and within-block taxa ordered by descending total abundance), the same
ordering convention used by the barplot, so a barplot and bubbleplot of the
same data present taxa in a consistent order.

## Styling

No title -- axis tick labels, a dot-size legend (relative abundance %), and
(if `group_rank` set) a color legend for the group.

## Example

```python
from microplotpy import load
from microplotpy.plots import taxa_bubbleplot

data = load("feature_table.tsv", "taxonomy.tsv", "metadata.tsv")
fig = taxa_bubbleplot(data, rank="Genus", group_rank="Phylum", top_n=25)
fig.savefig("bubbleplot.png", dpi=300, bbox_inches="tight")
```

## Literature

- Bubble/dot plot convention: `phyloseq::psmelt()` + ggplot bubble, MicrobiomeAnalyst dot plots.
- Complete-grid (zero-size, not omitted) convention for prevalence reading.
