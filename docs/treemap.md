# Taxa / functional treemap

`mp_taxa_treemap()` / `mp_function_treemap()` (R) -- `taxa_treemap()` /
`function_treemap()` (Python) -- hierarchical composition as nested
rectangles, a static publication-ready alternative to an interactive Krona
chart.

## Required input

Same as [mp_taxa_barplot()]/`taxa_barplot()` (taxonomic) or
[mp_function_barplot()]/`function_barplot()` (functional) -- reuses the
same `rank`/`group_rank`/filtering/tax-fix conventions as those, just
rendered as a treemap instead of a barplot.

## How it works

Area = relative abundance, **mean across samples** (not summed) so total
area is ~100% regardless of sample count -- comparable to the barplot's
per-sample relative abundance (%), just aggregated across the whole
dataset rather than shown per sample. Rectangles are grouped by
`group_rank` with a bordered outer region per group (bold label
top-left) and squarified-packed `rank` rectangles inside, colored by group.
Set `group_rank = NULL`/`None` for a flat treemap (one color per taxon, no
grouping borders).

Every label is clipped to its own rectangle, so long names in narrow
rectangles are cropped rather than bleeding into neighboring cells or off
the figure -- there is no external legend; every rectangle is labeled
directly, the way Krona charts are read.

## Key parameters

| Param | Default | Effect |
|---|---|---|
| `rank` | `"Genus"` | Level shown as inner rectangles. |
| `group_rank` | `"Phylum"` | Upper level for grouping/coloring/bordering. `NULL`/`None` for flat coloring. |
| `top_n` | `NULL`/`None` | Individual taxa shown; rest pooled into `other_label`. `NULL` shows everything passing the filters. |
| `min_rel_abund`, `min_prevalence`, `detection` | 0 | Same filters as the barplot -- see `barplot.md`. |
| `fix_taxonomy` | `TRUE`/`True` | Same tax-fix fallback as the barplot. |

## Example

```python
from microplotpy import load
from microplotpy.plots import taxa_treemap, function_treemap

data = load("feature_table.tsv", "taxonomy.tsv", "metadata.tsv")
fig = taxa_treemap(data, rank="Genus", group_rank="Phylum")
fig.savefig("treemap.png", dpi=300, bbox_inches="tight")

# functional data (eggNOG-mapper shape)
fig2 = function_treemap(gene_counts, annotation, rank="KEGG_ko", group_rank="COG_category")
```

## Literature

- Krona: interactive hierarchical taxonomic/functional visualization (Ondov BD et al. 2011, *BMC Bioinformatics*, 12:385) -- this treemap targets the same use case (readable hierarchical composition in one view) as a static figure.
- Treemap layout: Bruls M, Huizing K, van Wijk J. 2000, "Squarified Treemaps", *Eurographics*.
- R implementation: `treemapify` (Wilkins D). Python: `squarify`.
