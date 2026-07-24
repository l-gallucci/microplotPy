# Beta diversity

`mp_beta_diversity_plot()` (R) -- `beta_diversity_plot()` (Python) --
Bray-Curtis or Jaccard dissimilarity, PCoA or NMDS ordination, colored by a
grouping variable with 95% confidence ellipses, PERMANOVA annotated on the
plot. The standard beta-diversity ordination figure in amplicon/marker-gene
studies.

## Required input

Standard tidy input. `metadata.tsv` must include the categorical grouping
column (`group_var`) used for coloring, ellipses, and PERMANOVA.

## Key parameters

| Param | Default | Effect |
|---|---|---|
| `method` | `"bray"` | `"bray"` (Bray-Curtis, abundance-weighted) or `"jaccard"` (presence/absence only -- computed on the binary/turnover pattern, per ecological convention: this is *why* Jaccard and Bray-Curtis are usually reported together, one weighted by abundance and one not). |
| `ordination` | `"pcoa"` | `"pcoa"` (classical/metric MDS, axes = % variance explained) or `"nmds"` (non-metric MDS, annotated with stress). |
| `show_ellipse` | `TRUE`/`True` | 95% confidence ellipse per group. |
| `permanova` | `TRUE`/`True` | Run and annotate a PERMANOVA test for group separation. |
| `permutations` | 999 | Permutations for the PERMANOVA p-value. |

No phylogenetic (UniFrac) distance -- Bray-Curtis/Jaccard only, consistent
with the rest of this library (tree-based analysis is treated as a separate
step).

## PERMANOVA

R uses `vegan::adonis2()` directly. Python implements the same one-way
pseudo-F permutation test from first principles (Anderson 2001): observed
pseudo-F from the distance-based sums of squares, compared against pseudo-F
from many random relabelings of the group variable (distance matrix held
fixed) to get a p-value -- no extra dependency (`scikit-bio` avoided
deliberately; this is the same core math `adonis2`/`scikit-bio`'s
`permanova` use).

## Ordination math

PCoA is classical (metric) MDS: double-centre the squared distance matrix,
eigendecompose, scale eigenvectors by `sqrt(eigenvalue)` -- R uses
`stats::cmdscale()`, Python implements the identical linear-algebra
directly (no extra dependency). NMDS is non-metric MDS (Kruskal's stress
minimization) -- R uses `vegan::metaMDS()`, Python uses
`sklearn.manifold.MDS(metric_mds=False)`.

## Styling

No plot title in the ordinary sense -- PERMANOVA R² and p-value (and NMDS
stress) are essential statistical content here, not decorative chrome, so
they're shown as a subtitle/title-left annotation even though every other
plot in this library omits titles.

## Example

```python
from microplotpy import load
from microplotpy.plots import beta_diversity_plot

data = load("feature_table.tsv", "taxonomy.tsv", "metadata.tsv")
fig1 = beta_diversity_plot(data, group_var="Group", method="bray", ordination="pcoa")
fig1.savefig("beta_pcoa.png", dpi=300, bbox_inches="tight")

fig2 = beta_diversity_plot(data, group_var="Group", method="bray", ordination="nmds")
fig2.savefig("beta_nmds.png", dpi=300, bbox_inches="tight")
```

## Literature

- PERMANOVA: Anderson MJ. 2001, "A new method for non-parametric multivariate analysis of variance", *Austral Ecology*, 26(1):32-46.
- NMDS: Kruskal JB. 1964, "Nonmetric multidimensional scaling", *Psychometrika*, 29(1):1-27.
- Bray-Curtis: Bray JR, Curtis JT. 1957, *Ecological Monographs*, 27(4):325-349.
