# PRA-B---Localization-paper

## Code for the PRA paper

This repository contains the MATLAB code used to study **continuous-time quantum walks (CTQWs)** on three families of graphs and to compare a numerically computed quantity, the **dynamical inverse participation ratio (IPR)**, with closed-form analytical formulas.

The main script is:

```matlab
IPR_limiting_combined.m
```

The script produces the two figures described below and is fully self-contained: no external datasets or `.mat`/`.csv` files are required.

## Requirements

- MATLAB
- No external datasets are required.
- No additional data files are required.

The graphs and all numerical data used by the calculations are generated programmatically by the MATLAB script.

## Running the code

1. Clone or download this repository.
2. Open MATLAB.
3. Set the MATLAB current directory to the repository directory.
4. Run:

```matlab
IPR_limiting_combined
```

The script generates the two figures described below.

## 1. The three graph families

All three graphs are built from `n`-vertex complete graphs (cliques, `K_n`) glued together in different ways:

- **Barbell graph `B_n`** (`build_barbell`) — two cliques `K_n`, joined by a single bridge edge between one vertex of each clique. Total vertices: `2n`.

- **Star of Cliques, Variant 1 — "full connection"** (`build_star_v1`) — one hub vertex connected to every vertex of `n` separate `K_n` cliques. Total vertices: `1 + n²`.

- **Star of Cliques, Variant 2 — "single connection"** (`build_star_v2`) — the same layout, but the hub connects to only one vertex per clique (that vertex acts as a bridge vertex for the clique). Total vertices: `1 + n²`.

## 2. Core computation pipeline

For each adjacency matrix `A`, the script:

1. **Normalizes the adjacency matrix** (`normalize_adj`):
   `M = D^{-1/2} A D^{-1/2}`, i.e. a symmetric, degree-normalized version of the adjacency matrix.

2. **Computes the limiting (Cesàro-time-averaged) distribution** (`limiting_distribution`): diagonalizes `M`, groups eigenvectors by (numerically) repeated eigenvalues into eigenspaces, builds the orthogonal projector `P_E` onto each eigenspace, and accumulates

   `Π = Σ_E (P_E)∘(P_E)`

   where `∘` denotes the elementwise product. The entry `Π_ij` is the long-time-averaged probability of finding a quantum walker at vertex `i` given that it started at vertex `j`.

3. **Computes the dynamical IPR** (`dynamical_ipr`) for a chosen starting vertex `j`:

   `IPR_j = Σ_i Π_ij²`.

   This measures how localized (large IPR) versus how spread out (small IPR) the long-time walk distribution is when starting at vertex `j`.

## 3. Figure 1 — IPR versus graph size

For

```matlab
n = 4, 6, 8, ..., 70
```

the script computes the numerical dynamical IPR at specific vertices of interest (clique-interior vertices, bridge vertices, and hub vertices) for all three graph families.

These numerical values are compared with analytical/exact formulas derived from the known spectra of the graphs, including formulas such as

```text
1 - 4/n + 1/n²
```

for a clique-interior vertex,

```text
(n⁴ + 2n² + 5)/(n+1)⁴
```

for the Variant-1 hub, and

```text
1/4
```

for the Variant-2 hub.

**Figure 1** is a three-panel plot, one panel for each graph family. Solid lines represent the theoretical formulas and markers represent numerical values obtained by direct diagonalization. The figure therefore provides a numerical check of the analytical predictions.

## 4. Figure 2 — Full limiting-distribution heatmaps

For a fixed size

```matlab
n = 6
```

the script computes the full matrix `Π` for each of the three graph families and displays it as a heatmap using `imagesc`.

**Figure 2** is a three-panel plot showing, for every pair of vertices `(i,j)`, the long-time probability of transition from `j` to `i`.

The heatmaps visualize the block/community structure of the graphs: cliques appear as bright blocks, while bridge and hub connections produce fainter off-block regions.

## 5. Helper functions

The main script contains the following functions:

| Function | Purpose |
|---|---|
| `build_barbell(n)` | Builds the adjacency matrix for the barbell graph |
| `build_star_v1(n)` | Builds the adjacency matrix for Star-of-Cliques V1 |
| `build_star_v2(n)` | Builds the adjacency matrix for Star-of-Cliques V2 |
| `normalize_adj(A)` | Performs symmetric degree-normalization of the adjacency matrix |
| `limiting_distribution(M)` | Computes the Cesàro-limit transition-probability matrix `Π` via spectral projectors |
| `dynamical_ipr(Pi, j)` | Computes the dynamical IPR for a walk started at vertex `j` |
| `clamp01(x)` | Utility that clips a value into `[0,1]`; defined but not called in the current script |

## 6. Data and reproducibility

No external data are used or required.

All graphs are generated programmatically inside `IPR_limiting_combined.m` using the combinatorial definitions of the three graph families. The numerical parameters used in the two parts of the script are:

- `nvals = 4:2:70` for Figure 1.
- `n_fig = 6` for Figure 2.

The script does not read external datasets or `.mat`/`.csv` files and does not save output files automatically. The figures are displayed as MATLAB figure windows when the script is run.

## 7. Citation

If you use this code in your research, please cite the accompanying paper:

> **[Paper citation to be added when the final publication details are available.]**

The repository can be updated with the final bibliographic information and DOI when available.

## 8. License

This code is distributed under the MIT License. See the [`LICENSE`](LICENSE) file for details.

## 9. Repository contents

At present, the repository contains the MATLAB script

```text
IPR_limiting_combined.m
```

together with this README and the MIT license.
