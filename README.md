# microplotpy

Literature-grounded, publication-ready microbial ecology plots from plain
tidy input (feature table + taxonomy + metadata), with an upload-time
validator that reports exactly what is wrong with your data, and a
Shiny-for-Python app for interactive use. Python counterpart to the R
package [`microplotr`](https://github.com/l-gallucci/microplotR) — same
input format, same validator rules, same plot catalog, different plotting
engine (matplotlib/seaborn here vs ggplot2 on the R side).

See [`data-format.md`](data-format.md) for the required input files. Full
documentation (one page per plot: input spec, parameters, literature
grounding) is an [mkdocs-material](https://squidfunk.github.io/mkdocs-material/)
site built from `docs/` — run `mkdocs serve` for a live preview, or
`mkdocs build` to render static HTML into `site/` (e.g. for GitHub Pages
once this repo is published).

## Plot catalog

- **Taxonomy**: `taxa_barplot()` (nested-legend stacked barplot),
  `taxa_heatmap()` (CLR/log10, clustered), `taxa_bubbleplot()`,
  `taxa_gradient_plot()` (abundance over a continuous gradient),
  `alpha_diversity_plot()`, `beta_diversity_plot()` (PCoA/NMDS +
  PERMANOVA), `taxa_treemap()`.
- **Functional profiles** (eggNOG-mapper/KofamScan-shaped):
  `function_barplot()`, `function_heatmap()`, `function_treemap()` — thin
  wrappers around the taxa engine.
- **Metagenomics**: `mag_quality_plot()`/`mag_quality_distribution()`
  (CheckM/CheckM2), `assembly_nx_plot()`/`assembly_summary_barplot()`
  (QUAST-shaped).

Every plot resolves missing/unclassified annotation via `tax_fix()`,
supports abundance/prevalence filtering via `filter_taxa()`, and returns a
plain `matplotlib.Figure` you can keep mutating (`fig.axes[0].set_title(...)`
etc.) — nothing here is a static image.

## Install (dev)

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[plots,diversity,app,dev,docs]"
pytest
```

## Shiny app

```bash
shiny run app/app.py
```

Four tabs (Taxonomy, Functional profile, MAG quality, Assembly QC), each:
upload the required file(s) → validator shows errors/warnings inline →
pick a plot type and parameters → view/download the plot as PNG or SVG at
a width/height you choose.

## Quick usage (library)

```python
from microplotpy import load, validate
from microplotpy.plots import taxa_barplot

data = load(
    "examples/example_valid/feature_table.tsv",
    "examples/example_valid/taxonomy.tsv",
    "examples/example_valid/metadata.tsv",
)
report = validate(data, gradient_column="Depth_m", group_column="Group")
print(report.is_valid, report.errors, report.warnings)

fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=10)
fig.savefig("barplot.png", dpi=300, bbox_inches="tight")
```

Try validation against `examples/example_broken_*` to see what it catches
(ID mismatches, negative counts, missing taxonomy columns, duplicate
feature IDs).
