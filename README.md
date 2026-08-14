# What This Code Does

This MATLAB script studies **continuous-time quantum walks (CTQWs)** on three
families of graphs, and compares a *numerically computed* quantity called the
**dynamical inverse participation ratio (IPR)** against *closed-form analytical
formulas*. It produces two figures.

## 1. The three graph families

All three graphs are built from `n`-vertex complete graphs (cliques, `K_n`)
glued together in different ways:

- **Barbell graph `B_n`** (`build_barbell`) — two cliques `K_n`, joined by a
  single bridge edge between one vertex of each clique. Total vertices: `2n`.
- **Star of Cliques, Variant 1 — "full connection"** (`build_star_v1`) — one
  hub vertex connected to *every* vertex of `n` separate `K_n` cliques.
  Total vertices: `1 + n²`.
- **Star of Cliques, Variant 2 — "single connection"** (`build_star_v2`) —
  same layout, but the hub connects to only *one* vertex per clique (that
  vertex acts as a "bridge" vertex for the clique). Total vertices: `1 + n²`.

## 2. Core computation pipeline (for each graph)

For each adjacency matrix `A`, the script:

1. **Normalizes** the adjacency matrix (`normalize_adj`):
   `M = D^{-1/2} A D^{-1/2}`, i.e. a symmetric, degree-normalized version of
   the adjacency matrix (the natural Hamiltonian/generator for a CTQW on a
   graph).
2. **Computes the limiting (Cesàro-time-averaged) distribution**
   (`limiting_distribution`): diagonalizes `M`, groups eigenvectors by
   (numerically) repeated eigenvalues into eigenspaces, builds the orthogonal
   projector `P_E` onto each eigenspace, and accumulates
   `Π = Σ_E (P_E)∘(P_E)` (elementwise square of each projector, summed over
   distinct eigenvalues `E`). The entry `Π_ij` is the long-time-averaged
   probability of finding a quantum walker at vertex `i` given it started at
   vertex `j`.
3. **Computes the dynamical IPR** (`dynamical_ipr`) for a chosen starting
   vertex `j`: `IPR_j = Σ_i Π_ij²`. This measures how *localized* (large IPR,
   walker stays put) vs. how *spread out* (small IPR, walker delocalizes
   over many vertices) the long-time walk distribution is when starting at
   vertex `j`. It is the quantum-walk analogue of the IPR used to
   characterize localization in condensed matter physics.

## 3. Part 1 — IPR vs. graph size `n`

For `n = 4, 6, 8, ..., 70`, the script computes the **numerical** dynamical
IPR at specific vertices of interest (clique-interior vertices, bridge
vertices, hub vertices) for all three graph families, and compares them to
**analytical/exact formulas** derived from the graphs' known spectra (e.g.
`1 - 4/n + 1/n²` for a clique-interior vertex, `(n⁴+2n²+5)/(n+1)⁴` for the
Variant-1 hub, `1/4` for the Variant-2 hub, etc.).

The result is **Figure 1**, a 3-panel plot (one panel per graph family)
showing solid lines (theoretical formulas) overlaid with markers (numerical
values) as functions of `n` — a visual check that the analytical predictions
match direct numerical diagonalization.

## 4. Part 2 — Full limiting-distribution heatmaps

For a fixed size (`n = 6`), the script computes the full matrix `Π` for each
of the three graphs and displays it as a **heatmap** (`imagesc`) — **Figure
2**, a 3-panel plot showing, for every pair of vertices `(i, j)`, the
long-time probability of transition from `j` to `i`. This visualizes the
block/community structure of each graph (cliques appear as bright blocks,
bridges/hub connections as faint off-block regions).

## 5. Helper functions summary

| Function | Purpose |
|---|---|
| `build_barbell(n)` | Builds adjacency matrix for the barbell graph |
| `build_star_v1(n)` | Builds adjacency matrix for Star-of-Cliques V1 (hub fully connected to each clique) |
| `build_star_v2(n)` | Builds adjacency matrix for Star-of-Cliques V2 (hub connected to one vertex per clique) |
| `normalize_adj(A)` | Symmetric degree-normalization of adjacency matrix |
| `limiting_distribution(M)` | Computes the exact Cesàro-limit transition-probability matrix `Π` via spectral projectors |
| `dynamical_ipr(Pi, j)` | Computes the dynamical IPR for a walk started at vertex `j` |
| `clamp01(x)` | Utility to clip a value into `[0, 1]` (defined but not actually called anywhere in the script) |

---

## What about data?

**No external data is used or required.** This script is fully
**self-contained**:

- It does not read any files, datasets, `.mat` files, or `.csv` files.
- All "data" — the graphs (adjacency matrices) — are **generated
  programmatically** inside the script itself, using pure combinatorial rules
  (`ones(n) - eye(n)` for a complete graph block, plus explicit bridge/hub
  edges).
- The only numerical inputs are the graph-size parameters `nvals = 4:2:70`
  (for Part 1) and `n_fig = 6` (for Part 2), both hard-coded at the top of
  each section.
- The only outputs are the two MATLAB figures (`figure(...)` windows) — the
  script does not save any files, write results to disk, or export data,
  unless you add that yourself (e.g. `saveas(gcf, 'fig1.png')`).

So there's nothing to upload or supply — you can run the script as-is and it
will regenerate everything from scratch.
# PRA-B---Localization-paper
