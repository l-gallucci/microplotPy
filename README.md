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
site built from `docs/`, served via GitHub Pages once enabled for this repo.

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

## Install

```bash
pip install "microplotpy[plots,diversity] @ git+https://github.com/l-gallucci/microplotPy.git"
```

Drop `[plots,diversity]` for the bare library (validator only, no plotting
deps). Not yet on PyPI — once published there, this becomes a plain
`pip install microplotpy[plots,diversity]`.

## Quick usage

Point the loader at your own tidy files (see
[`data-format.md`](data-format.md) for the required columns):

```python
from microplotpy import load, validate
from microplotpy.plots import taxa_barplot

data = load("feature_table.tsv", "taxonomy.tsv", "metadata.tsv")
report = validate(data, gradient_column="Depth_m", group_column="Group")
print(report.is_valid, report.errors, report.warnings)

fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=10)
fig.savefig("barplot.png", dpi=300, bbox_inches="tight")
```

## Shiny app

The app (`app/app.py`) isn't packaged for `pip install` yet — it needs the
repo cloned:

```bash
git clone https://github.com/l-gallucci/microplotPy.git && cd microplotPy
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[plots,diversity,app]"
shiny run app/app.py
```

Four tabs (Taxonomy, Functional profile, MAG quality, Assembly QC), each:
upload the required file(s) → validator shows errors/warnings inline →
pick a plot type and parameters → view/download the plot as PNG or SVG at
a width/height you choose.

## Contributing / development

```bash
git clone https://github.com/l-gallucci/microplotPy.git && cd microplotPy
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[plots,diversity,app,dev,docs]"
pytest
```

`examples/` ships several bundled datasets, including `example_valid/` (used
in the docs) and `example_broken_*/` sets that each trip a different
validator check (ID mismatches, negative counts, missing taxonomy columns,
duplicate feature IDs):

```python
data = load(
    "examples/example_valid/feature_table.tsv",
    "examples/example_valid/taxonomy.tsv",
    "examples/example_valid/metadata.tsv",
)
```

To rebuild the documentation site after a change: `mkdocs build`.
