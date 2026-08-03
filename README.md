# Trophic Maelstrom 
## A High Performance Simulation and Data-Analysis Pipeline For Desertification Transitions. ⚡⚡

A scientific HPC SPDE (Stochastic Partial Differential Equation) codebase for spatial ecological simulations: the Rietkerk vegetation-water dryland model (extended to a 3-species vegetation-grazer-predator trophic chain) and directed-percolation
(DP) model (implemented as the Reggeon Field Theoretic Equation, similarly extended to 1/2/3 species trophic chains) as used as simpler analogues for studying critical phenomena. The codebase spans three layers:
C++/CUDA simulation engines, a Python/Bash layer that collates simulation data, and then reorganises and statistically post-processes their raw output, and a Python visualisation layer that turns the result into heatmaps, videos, and summary plots.

## Architecture

![Pipeline architecture flowchart](architecture_flowchart.png)

Data flows top to bottom through three layers, each writing into a predictable directory tree that the next
layer consumes. Two feedback loops run the other way: `Utilities/synthetic_clump_generator.py` seeds
synthetic initial-condition frames back into `Input/` for the C++ layer to read, and
`Utilities/gmm_classifier.py` is called live, in-process, from `Percolation_FD`'s C++ code (under a
`-DCRITGMM` build) via an embedded CPython bridge rather than through a file on disk.

### 1. C++ / CUDA simulation layer

Three model implementations, each independently compiled (no Makefile -- see each README for the exact
`g++`/`nvcc` invocations and compile-time macros):

- **[`simulations/Rietkerk_FD/`](simulations/Rietkerk_FD/README.md)** -- finite-difference stochastic
  implementation of the Rietkerk vegetation-water model (1/2/3-species trophic chain) (CPU-validated ✅, CUDA accelerated path contains **known bugs** ❌ ).
- **[`simulations/Rietkerk_FPE/`](simulations/Rietkerk_FPE/README.md)** -- Fokker-Planck-equation variant of
  the same model, adding an advection-diffusion step (GPU-validated ✅; CPU-only path is more experimental and may contain unverified bugs ⚠️⚠️).
- **[`simulations/Percolation_FD/`](simulations/Percolation_FD/README.md)** -- multi-species
  scale free model, implemented as the Reggeon Field Theoretic Equation, and belonging to the Directed Percolation (DP) universality class (1/2/3-species trophic chain) (CPU-validated ✅, CUDA accelerated path contains **known bugs** ❌ ).

All three write raw per-replicate CSV output into an ad hoc directory tree under `simulations/Data/`.

### 2. Python/Bash reorganisation & analysis layer

**[`simulations/Utilities/`](simulations/Utilities/README.md)** turns that raw C++ output across different devices (requires rsync) into a predictable, analysis-ready directory structure and computes a battery of spatial/temporal statistics over it, all
CPU/GPU-agnostic ✅ ([`slick.py`](https://github.com/HAL9K000/slick) routes numerics through `cupy` when available, `numpy` otherwise). Functionality
provided:

- **Reorganisation**: copies and renames raw per-replicate CSV files (`reorganise_dir.py` for Frame
  snapshots, `reorganise_prelims_dir.py` for Prelim time series) into a predictable
  `{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/...` tree, with collision-safe replicate handling.
- **Per-snapshot/per-directory statistics**: replicate mean/std and survival counts, per-cell replicate
  means, 2D spatial FFT power spectra.
- **Spatial pattern analysis**: clustering (KMeans/GMM/threshold), cluster-size and macro cluster
  statistics, KDE-based potential-well/bistability analysis.
- **Spatio-temporal synchrony analysis**: 2D FFT cross-correlation (NCC/ZNCC), mutual information (MI/AMI),
  and bivariate Moran's I between snapshots at different timepoints, plus harmonic-frequency extraction from
  the resulting correlation curves -- GPU-accelerated for large grids.
- **Batch orchestration & remote data collation**: `dreamcatcher_automata.bash`/`prelims_automata.bash` run
  the above reorganisation and post-processing across many parameter combinations unattended, and (via
  `expect_commands.sh`) can first `rsync` raw output in from multiple remote machines/HPC clusters into a
  single host before reorganising it.
- **Synthetic initial conditions**: `synthetic_clump_generator.py` generates synthetic hexagonal-lattice
  seed frames for the C++ simulators.
- **Live C++ bridge**: `gmm_classifier.py`, embedded directly into `Percolation_FD`'s C++ process for
  real-time GMM cluster-statistics computation during a running simulation.

### 3. Data visualisation layer

**[`simulations/Data_Processing/`](simulations/Data_Processing/README.md)** consumes the reorganised data and
derived statistics from the layer above and renders them. Functionality provided (via `frame_builder.py`):

- **Frame (spatial snapshot) visualisation**: per-species heatmap PNGs and stitched-together videos, in
  several groupings (all replicates combined, one video per replicate, across-`a` at fixed `T`, with a GAMMA
  interaction-field overlay for spreading-test runs).
- **Frame-derived summary plots**: equilibrium state vs. control parameter `a`, FFT power-spectrum heatmaps,
  potential-well/bistability plots, cluster-size distributions and macro cluster-statistic time series.
- **Prelims (time-series) summary plots**: time-series plots with optional power-law/decay fit overlays,
  equilibrium-state summary and violin plots, phase-space trajectory plots, and cross-prefix/parameter
  comparison plots.

A companion Jupyter notebook, `Data_Processing/frame_builder_walkthrough.ipynb`, provides a full runnable,
documented walkthrough of every one of these functions -- see the `Data_Processing/README.md` for details.

## Repository layout at a glance

```
simulations/
  Rietkerk_FD/        C++/CUDA -- vegetation-water model
  Rietkerk_FPE/        C++/CUDA -- Fokker-Planck (advection-diffusion) variant
  Percolation_FD/      C++/CUDA -- directed-percolation model
  Utilities/            Python -- reorganisation, statistics, batch orchestration
  Data_Processing/       Python -- visualisation (frame_builder.py + notebook)
  Data/                 raw + reorganised simulation output
  Input/                 synthetic/seed initial-condition frames
Images/                 generated heatmaps, videos, summary plots
```

There is no top-level build system or package manager ⚡⚡ (as of yet) -- each layer is built/run independently; see the
README linked for each component above for exact commands.
