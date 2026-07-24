# microplotpy

Literature-grounded, publication-ready microbial ecology plots from plain
tidy input (feature table + taxonomy + metadata), with an upload-time
validator that reports exactly what is wrong with your data, and a
Shiny-for-Python app for interactive use.

Python counterpart to the R package `microplotr` (sibling repo) — same
input format, same validator rules, same plot catalog, different plotting
engine (matplotlib/seaborn here vs ggplot2 on the R side).

## Install

```bash
pip install -e ".[plots,diversity,app,dev]"
```

## Quick start

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

See [Data format](data-format.md) for the required input files (Feature_ID
/ Sample_ID conventions, validator rules), and the plot pages in the nav
for input specs, parameters, and the literature each plot is grounded in.

## Shiny app

```bash
shiny run app/app.py
```

Four tabs (Taxonomy, Functional profile, MAG quality, Assembly QC), each:
upload the required file(s) → validator shows errors/warnings inline →
pick a plot type and parameters → view/download the plot as PNG or SVG at
a size you choose.

## Plot catalog

**Taxonomy** — [Barplot](barplot.md), [Heatmap](heatmap.md),
[Bubbleplot](bubbleplot.md), [Gradient](gradient.md),
[Alpha diversity](alpha_diversity.md), [Beta diversity](beta_diversity.md),
[Treemap](treemap.md).

**Functional profiles** (eggNOG-mapper/KofamScan-shaped) —
[Function profile](function_profile.md) (barplot/heatmap/treemap wrappers
around the taxa engine).

**Metagenomics** — [MAG quality](mag_quality.md) (CheckM/CheckM2),
[Assembly QC](assembly.md) (QUAST-shaped Nx curve + summary).

Every plot resolves missing/unclassified annotation via `tax_fix()`,
supports abundance/prevalence filtering via `filter_taxa()`, and returns a
plain `matplotlib.Figure` you can keep mutating — nothing here is a static
image.
