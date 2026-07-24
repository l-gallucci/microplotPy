# Taxa heatmap

`mp_taxa_heatmap()` (R) / `taxa_heatmap()` (Python) -- clustered heatmap of
the top-N taxa at a chosen rank x samples, CLR- or log10-transformed.

## Required input

Same tidy input as the barplot. Only `feature_table.tsv` + `taxonomy.tsv` at
the chosen `rank` are used (metadata not required unless you filter first).

## Key parameters

| Param | Default | Effect |
|---|---|---|
| `rank` | `"Genus"` | Taxonomy level shown as rows. |
| `top_n` | 25 | Taxa shown (ranked by total relative abundance among taxa passing the filters below). `NULL`/`None` shows every taxon that passes the filters, uncapped. |
| `min_rel_abund` | 0 | Minimum mean relative abundance (%) across all samples for a taxon to be shown. `0` disables. |
| `min_prevalence` | 0 | Minimum fraction of samples (`0-1`), or an absolute count if `> 1`, in which a taxon must be detected to be shown. `0` disables. |
| `detection` | 0 | Relative-abundance (%) threshold above which a taxon counts as "detected" in a sample, for `min_prevalence`. |
| `transform` | `"clr"` | `"clr"` (centered log-ratio, per sample) or `"log10"` (`log10(x + pseudocount)`). |
| `cluster_rows`, `cluster_cols` | `TRUE`/`True` | Hierarchical clustering (UPGMA / average linkage) on each axis. |
| `show_dendrogram` | `TRUE`/`True` | Draw dendrograms on both margins. |
| `fix_taxonomy` | `TRUE`/`True` | Same tax-fix step as the barplot -- see `barplot.md`. |

## Filtering

`min_rel_abund`/`min_prevalence` narrow the candidate taxa first (a taxon
must pass both if both are set), then `top_n` picks the most abundant among
survivors -- see `barplot.md` for the full rationale. Unlike the barplot,
taxa that don't make the cut are simply not shown (no "Other" row).

## Distance choice

Samples are clustered by **Bray-Curtis dissimilarity on relative
abundance** -- the community-ecology standard for compositional microbiome
data (used by default in `vegan`, `ampvis2::amp_heatmap`, QIIME2 beta
diversity). Taxa are clustered by **Euclidean distance on the transformed
(CLR or log10) values actually shown**, matching the `pheatmap` /
`ComplexHeatmap` default for the displayed matrix.

## Transform choice

CLR (centered log-ratio) is the compositional-data-appropriate transform for
microbiome relative abundances (Gloor et al. 2017, *Frontiers in
Microbiology* -- "microbiome datasets are compositional") and is the
default. `log10(x + pseudocount)` is offered as a simpler, more familiar
alternative for exploratory use.

## Styling

No title -- axis tick labels (taxa/samples) and a single colorbar legend
only. Row/column order follows the clustering, not alphabetical/abundance
order, so the dendrogram is directly readable against the tile grid.

## Example

```python
from microplotpy import load
from microplotpy.plots import taxa_heatmap

data = load("feature_table.tsv", "taxonomy.tsv", "metadata.tsv")
fig = taxa_heatmap(data, rank="Genus", top_n=25, transform="clr")
fig.savefig("heatmap.png", dpi=300, bbox_inches="tight")
```

## Literature

- Bray-Curtis clustering convention: `ampvis2::amp_heatmap`, `vegan::vegdist`.
- CLR transform for compositional microbiome data: Gloor GB et al. 2017, *Frontiers in Microbiology*, 8:2224.
- General clustered-heatmap layout: `pheatmap`, `ComplexHeatmap` (Gu Z et al. 2016, *Bioinformatics*).
