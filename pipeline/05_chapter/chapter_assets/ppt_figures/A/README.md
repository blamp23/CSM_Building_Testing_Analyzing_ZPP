# Figure A — Where the models differ

Two outputs, both built by `pipeline/05_chapter/scripts/11_master_figures.py`.

## Chapter figure

**`A_where_models_differ.{pdf,png}`** — two panels, 16:9.

- **A1.** Stacked bars of notable representatives (|centered log₂FC| ≥ 1,
  carried by both cells, loop-flagged excluded) by metabolic class. Left
  group: five F1 bars (one per hpf stage), pooled over the three cross-regime
  contrasts at that stage. Right group: four F2 bars (one per adjacent-stage
  transition), pooled over the three regimes. Class colours match SA.
- **A2.** Centered PCA of the 15 cells (PC1 vs PC2). Regime by marker, hpf
  by colour saturation, 120 hpf cells bold. Values from
  `../../comparison_tables/pca_cells_centered.csv`.

## Supplementary master

**`SA_contrast_class_bars.{pdf,png}`** — 27 panels, one per individual
contrast (15 F1 + 12 F2), sorted by hpf. Same class colours as the chapter
view. Every pooled bar in A1 is the sum of three of these panels.

## Supplementary table

**`class_map.csv`** — the subsystem → class assignment used by both figures.
123 rows, one per subsystem.

**`SA_counts.csv`** — the raw count-by-class-per-contrast table that A1
sums over.

## Caveat

A1 is pooled over three contrasts per bar. SA shows the same signal
holds (or doesn't) in each contrast individually — a reader asking "is
lipid dominance at 96 hpf a BL-vs-D artifact?" should check SA, where
the same ~lipid height is visible in all three 96 hpf contrasts.
