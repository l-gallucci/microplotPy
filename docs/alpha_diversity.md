# Alpha diversity

`mp_alpha_diversity_plot()` (R) -- `alpha_diversity_plot()` (Python) --
Observed richness, Shannon, and Simpson diversity per sample, boxplot +
jitter by a grouping variable, one panel per metric, with a statistical
test annotated. The standard alpha-diversity comparison figure in
amplicon/marker-gene studies.

## Required input

Standard tidy input. `metadata.tsv` must include the categorical grouping
column you compare (`group_var`, ≥2 levels) -- checked by the validator's
`group_column` option (see `data-format.md`).

## Metrics

- **Observed**: richness -- count of taxa with nonzero abundance in a
  sample (R: `vegan::specnumber()`; Python: count of nonzero features).
- **Shannon**: `-sum(p * log(p))` over relative proportions `p` (natural
  log, matching `vegan::diversity(index = "shannon")`).
- **Simpson**: Gini-Simpson index `1 - sum(p^2)` (matching
  `vegan::diversity(index = "simpson")`).

Faith's phylogenetic diversity is intentionally excluded -- it needs a
phylogenetic tree, out of scope for this tidy-flat-file library (see
`data-format.md`).

## Statistical test

`test = "auto"` (default) picks Wilcoxon/Mann-Whitney U for exactly 2
groups, Kruskal-Wallis for more than 2 -- both non-parametric, appropriate
for diversity indices which are rarely normally distributed. Set
`test = NULL`/`None` to disable. If a metric is literally constant across
every sample (zero variance -- e.g. identical richness), the p-value
annotation is silently omitted for that panel rather than showing a
meaningless result.

## Styling

No plot title -- bold per-metric facet labels, boxplot + jittered points
(jitter is x-only; the y-position always shows the real value, never
artificially perturbed), and the test p-value as a small annotation above
each panel.

## Example

```python
from microplotpy import load
from microplotpy.plots import alpha_diversity_plot

data = load("feature_table.tsv", "taxonomy.tsv", "metadata.tsv")
fig = alpha_diversity_plot(data, group_var="Group")
fig.savefig("alpha_diversity.png", dpi=300, bbox_inches="tight")
```

## Literature

- Shannon H. 1948, "A Mathematical Theory of Communication", *Bell System Technical Journal*.
- Simpson EH. 1949, "Measurement of Diversity", *Nature*, 163:688.
- Non-parametric group comparison for diversity indices: standard practice, e.g. QIIME2 `diversity alpha-group-significance`.
