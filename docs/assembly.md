# Assembly QC: Nx curve + summary barplot

`mp_assembly_nx_plot()` / `mp_assembly_summary_barplot()` (R) --
`assembly_nx_plot()` / `assembly_summary_barplot()` (Python) -- the standard
QUAST-style assembly quality-control figures.

## Required input

Two separate flat tables, not the 3-file taxonomy schema:

**`contig_lengths.tsv`** (for the Nx curve) -- one row per contig:

| Assembly_ID | Contig_ID | Length |
|---|---|---|
| Assembly_A | Assembly_A_contig_1 | 122604 |

**`assembly_summary.tsv`** (for the summary barplot) -- one row per
assembly, shaped like QUAST's `report.tsv`. QUAST's own column names are
accepted directly (normalized internally):

| Canonical column | QUAST name |
|---|---|
| `Assembly_ID` | `Assembly` |
| `N_contigs` | `# contigs` |
| `Largest_contig` | `Largest contig` |
| `Total_length` | `Total length` |
| `N50` | `N50` |
| `L50` | `L50` |
| `GC_percent` | `GC (%)` |

`Assembly_ID`, `N50`, `Total_length` are required (under either naming
convention); the rest are optional and only needed if you plot them via
`stat_col`.

## Validators

`mp_validate_contig_lengths()`/`validate_contig_lengths()`: required
columns; `Length` numeric and positive; no duplicate `Contig_ID` within an
`Assembly_ID`.

`mp_validate_assembly_summary()`/`validate_assembly_summary()`: required
columns; `N50`/`Total_length` numeric and positive; no duplicate
`Assembly_ID`.

## Nx curve

Cumulative contig-length distribution: contigs sorted longest-to-shortest,
x = cumulative % of total assembly length, y = contig length (log10 scale),
one step-line per assembly. Reading off x = 50 gives the assembly's N50 --
this is the standard QUAST plot, more informative than the single N50 number
alone since it shows the whole length distribution, not just one point.

## Summary barplot

A single QUAST `report.tsv` statistic (default `N50`) plotted per assembly
as a plain bar -- set `stat_col` to `"Total_length"`, `"N_contigs"`, etc. to
plot a different one.

## Styling

No title -- axis titles, tick labels, and (Nx curve only) a color legend
per assembly.

## Example

```python
import pandas as pd
from microplotpy.plots import assembly_nx_plot, assembly_summary_barplot

contigs = pd.read_csv("contig_lengths.tsv", sep="\t")
summary_tbl = pd.read_csv("assembly_summary.tsv", sep="\t")

fig1 = assembly_nx_plot(contigs)
fig1.savefig("assembly_nx.png", dpi=300, bbox_inches="tight")

fig2 = assembly_summary_barplot(summary_tbl, stat_col="N50")
fig2.savefig("assembly_summary.png", dpi=300, bbox_inches="tight")
```

## Literature

- QUAST: Gurevich A et al. 2013, *Bioinformatics*, 29(8):1072-1075.
- Nx curve / N50 as an assembly contiguity metric: standard since the original genome assembly papers (e.g. the human genome draft, IHGSC 2001).
