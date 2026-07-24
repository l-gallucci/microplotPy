# MAG quality plot

`mp_mag_quality_plot()` / `mp_mag_quality_distribution()` (R) --
`mag_quality_plot()` / `mag_quality_distribution()` (Python) -- the standard
first figure (and its companion distribution plot) in any metagenome-
assembled genome (MAG) recovery paper.

## Required input

A single flat table, one row per MAG/bin -- **not** the 3-file taxonomy
schema used by the marker-gene plots. Shaped like CheckM2's
`quality_report.tsv` or CheckM (v1)'s `qa`/`bin_stats` table; column names
from either tool are accepted directly (normalized internally):

| Canonical column | CheckM2 name | CheckM (v1) name |
|---|---|---|
| `Name` | `Name` | `Bin Id` |
| `Completeness` | `Completeness` | `Completeness` |
| `Contamination` | `Contamination` | `Contamination` |
| `Genome_Size` | `Genome_Size` | `Genome size (bp)` |
| `Contig_N50` | `Contig_N50` | `N50 (contigs)` |
| `GC_Content` | `GC_Content` | `GC` |

`Completeness`, `Contamination`, `Genome_Size` are required (under either
naming convention); everything else is optional and only used if you pass it
as `size_col`/`color_col` (e.g. a taxonomy column from GTDB-Tk for coloring).

## Key parameters (`mp_mag_quality_plot()` / `mag_quality_plot()`)

| Param | Default | Effect |
|---|---|---|
| `size_col` | `"Genome_Size"` | Column mapped to dot size. `NULL`/`None` for uniform dots. |
| `color_col` | `NULL`/`None` | Column mapped to dot color (e.g. a taxon/group column). |
| `show_mimag_thresholds` | `TRUE`/`True` | Draw the MIMAG quality-tier reference lines/shading (see below). |

`mp_mag_quality_distribution()` / `mag_quality_distribution()` takes the
same table and a `bins` count; no `size_col`/`color_col` (it's a histogram,
not a scatter).

## MIMAG quality tiers

Reference lines (and a shaded high-quality region) follow the MIMAG standard
(Bowers et al. 2017, *Nature Biotechnology* -- "Minimum information about a
single amplified genome and a metagenome-assembled genome"):

- **High quality**: completeness ≥ 90%, contamination < 5%.
- **Medium quality**: completeness ≥ 50%, contamination < 10%.
- Below that: low quality.

## Styling

No title -- axis titles, tick labels, and (if `color_col` set) a color
legend, plus a size legend if `size_col` set. `mp_mag_quality_distribution()`
uses per-panel bold labels ("Completeness (%)" / "Contamination (%)") the
same way the barplot/bubbleplot use bold facet strip labels -- not a plot
title.

## Example

```python
import pandas as pd
from microplotpy.plots import mag_quality_plot, mag_quality_distribution

mag = pd.read_csv("mag_quality.tsv", sep="\t")
fig1 = mag_quality_plot(mag, size_col="Genome_Size", color_col="Phylum")
fig1.savefig("mag_quality.png", dpi=300, bbox_inches="tight")

fig2 = mag_quality_distribution(mag)
fig2.savefig("mag_quality_distribution.png", dpi=300, bbox_inches="tight")
```

## Literature

- MIMAG quality-tier standard: Bowers RM et al. 2017, *Nature Biotechnology*, 35(8):725-731.
- CheckM2: Chklovski A et al. 2023, *Nature Methods*.
- CheckM (v1): Parks DH et al. 2015, *Genome Research*, 25(7):1043-1055.
