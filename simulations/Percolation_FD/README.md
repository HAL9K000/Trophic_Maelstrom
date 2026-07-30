# Percolation_FD Documentation

Finite-difference (Dornic-scheme) stochastic implementation of a multi-species stochastic partial-directed-
percolation (SPDP) model, extended to a 1/2/3-species trophic chain (vegetation → grazer → predator). It is
used as a structurally simpler analogue of the Rietkerk model (see `../Rietkerk_FD`) for studying
absorbing-state phase transitions and critical phenomena, since it drops the vegetation-water coupling
while keeping the same trophic-interaction and stochastic-integration machinery. Simulations run on a 2D
periodic lattice, parallelised over OpenMP threads (one thread per replicate/control-parameter combination).
Species interaction terms ("gamma", the local grazing/predation coupling) are computed each timestep either
directly over a neighbourhood (O(N²)) or via a cached FFT convolution (O(N² log N)) using FFTW3, and initial
frames can be generated from closed-form expressions parsed at runtime with ExprTk. Each run scans a range
of the control parameter `a` and writes CSV frames, per-timestep aggregate ("Prelim") data, and summary
order-parameter statistics to disk.

This document covers the two core implementation files — `multiSPDP.cpp`/`MultiSPDP.h` and the
`order_*stocDP_unity.cpp`/`order_*stocDP_burnin.cpp` drivers that link against them — and how to compile and
run them. The directory also contains a separate, simpler single-species DP implementation
(`SPDP.h`/`stochasticSPDP.cpp` and its interactive analysis drivers `bifurc_prob.cpp`,
`finite_scaling_delta.cpp`, `spreading_exp_test.cpp`, etc.), a CUDA kernel, and older scratch code (`New/`);
these are not covered here — consult those files directly if you need them.

For the repository-wide picture (data-output conventions, how this fits with `Rietkerk_FD`, `Rietkerk_FPE`,
and the `Utilities`/`Data_Processing` post-processing scripts), see the top-level `simulations` README and
`/CLAUDE.md`.

## Dependencies

- `g++` with C++23 support (scheduling scripts use `g++-14`)
- OpenMP (`-fopenmp`)
- [FFTW3](http://www.fftw.org/) — FFT-accelerated interaction-kernel convolution
  (`#include <fftw3.h>` in `MultiSPDP.h`). Not vendored — build/install it yourself and make sure your
  compiler/linker can find its headers and `libfftw3`/`libfftw3_threads`.
- [ExprTk](https://github.com/ArashPartow/exprtk) (header-only) — symbolic-expression parsing used for frame
  initialization (`#include <exprtk.hpp>` in `MultiSPDP.h`). Not vendored — put it on your include path.
- [Armadillo](http://arma.sourceforge.net/) — linear algebra / GMM-clustering support, optional but recommended. Only pulled in when compiled with `-DARMA`, but the scheduling scripts link `-larmadillo` unconditionally, so if using these scripts, it should be installed regardless.
- Only needed if compiling with `-DCRITGMM` (see below): a Python 3 interpreter with NumPy development
  headers, and `../Utilities/gmm_classifier.py` importable at runtime.
- GNU `screen`, only if using the local parallel job launchers in `Scheduling Bashes/`.
- A SLURM environment for HPC use, only if using the cluster job-array script.

None of FFTW3, ExprTk, or Armadillo is vendored in this repo — install them yourself.

## Core files

- `MultiSPDP.h` / `multiSPDP.cpp` — the shared model implementation, compiled into every binary. Provides:
  grid types over plain nested `std::vector`; frame-initialization routines; the `FFTW3_CentralPlanner`
  struct that caches FFT plans/kernels for the interaction-neighbourhood convolution; the RK4/Dornic
  stochastic integrators for the 1/2/3-species variants; both direct and FFT-accelerated gamma-calculation
  routines; CSV/frame I/O; and (optionally, see `-DCRITGMM` below) the embedded-Python GMM-clustering hooks.
- `MultiSPDP_constants_{1,2,3}Sp.h` — pulled in by `MultiSPDP.h` based on the `-DSPB` value at compile time;
  defines species-count-specific constants (total species count `Sp`, CSV frame/prelim/GMM-cluster column
  headers).
- `Debug.h` — small third-party (Nadeau Software, CC-BY 3.0) helper for querying peak/current memory (RSS)
  usage, included by `multiSPDP.cpp` for diagnostics.
- `order_{1,2,3}stocDP_unity.cpp`, `order_{2,3}stocDP_burnin.cpp` — the driver files. Each contains only a
  `main()` and is compiled together with `multiSPDP.cpp` to produce one executable. `unity` drivers start
  from a fresh (random-speckle) initial condition; `burnin` drivers read the initial frame from a
  previously-generated CSV instead.

## Compiling

There is no Makefile/CMakeLists — every run configuration is a standalone `g++` invocation, compiling
`multiSPDP.cpp` together with **one** `order_*stocDP_*` driver file. Which species-count model and which
initial-condition routine gets built is controlled by **compile-time macros**, not runtime flags. This
mirrors the options exposed by `Scheduling Bashes/twilight_screening.bash`:

```
Usage: twilight_screening.bash <Path/to/init_file.txt> [Optional: <SpB=3> <Init_Types=0>]
  SpB:        number of biotic species (default 3)
  Init_Types: 0 = Homogeneous Initialisation, 1 = Random Speckles, 2 = Burn-In Frames (default 0)
```

| Macro | Values | Meaning |
|---|---|---|
| `-DSPB=<N>` | `1`, `2`, `3` | Number of biotic species: vegetation only, +grazer, +grazer+predator. Selects which `order_*stocDP_*` driver file to compile against, and which `MultiSPDP_constants_<N>Sp.h` gets pulled in. |
| `-DINIT=<N>` | `0`, `1`, `2` | Initial condition: `0` = homogeneous frame, `1` = random speckles ("unity"), `2` = burn-in frame read from a file on disk. `0`/`1` require an `order_*stocDP_unity.cpp` driver (which explicitly refuses to build/run with `-DINIT=2`); `2` requires the corresponding `order_*stocDP_burnin.cpp` driver. |
| `-DARMA` | (flag) | Enables the Armadillo-dependent code paths in `MultiSPDP.h`. Required together with `-DCRITGMM`. |
| `-DCRITGMM` | (flag) | Enables the embedded-Python GMM-clustering pipeline — see below. Optional; off by default. |

### 1-species (vegetation only)

```bash
g++-14 -O3 -march=native -DSPB=1 -DINIT=0 \
  multiSPDP.cpp order_1stocDP_unity.cpp \
  -larmadillo -L${HOME}/.local/lib -lfftw3_threads -lfftw3 -lm -fopenmp \
  -o dp_1sp.out -std=c++23
```

### 2-species (vegetation + grazer)

```bash
g++-14 -O3 -march=native -DSPB=2 -DINIT=0 \
  multiSPDP.cpp order_2stocDP_unity.cpp \
  -larmadillo -L${HOME}/.local/lib -lfftw3_threads -lfftw3 -lm -fopenmp \
  -o dp_2sp.out -std=c++23
```

### 3-species (vegetation + grazer + predator)

```bash
g++-14 -O3 -march=native -DSPB=3 -DINIT=0 \
  multiSPDP.cpp order_3stocDP_unity.cpp \
  -larmadillo -L${HOME}/.local/lib -lfftw3_threads -lfftw3 -lm -fopenmp \
  -o dp_3sp.out -std=c++23
```

### Burn-in initial condition (`-DINIT=2`)

Swap the `unity` driver for the matching `burnin` one (2- and 3-species only):

```bash
g++-14 -O3 -march=native -DSPB=3 -DINIT=2 \
  multiSPDP.cpp order_3stocDP_burnin.cpp \
  -larmadillo -L${HOME}/.local/lib -lfftw3_threads -lfftw3 -lm -fopenmp \
  -o dp_3sp_burnin.out -std=c++23
```

### CRITGMM (embedded-Python GMM clustering)

Setting `-DCRITGMM -DARMA` also compiles in `py_GMM_embedder.hpp`, which embeds a CPython interpreter,
imports `../Utilities/gmm_classifier.py`, and calls its `gen_clustered()` function live during frame
generation to compute Gaussian-mixture cluster statistics (mean cluster size, number of clusters, occupied
site fraction, etc.) for each species at every save point.

```bash
g++-14 -O3 -march=native -DSPB=3 -DINIT=0 -DCRITGMM -DARMA \
  multiSPDP.cpp order_3stocDP_unity.cpp \
  -L${HOME}/.local/lib -lfftw3_threads -lfftw3 -lm -fopenmp \
  $(python3-config --includes) -I$(python3 -c "import numpy; print(numpy.get_include())") \
  $(python3-config --ldflags --embed) \
  -std=c++23 -o dp_3sp_critgmm.out
```
This isn't wired into any of the Scheduling Bashes scripts as an active option — treat it as experimental, for the moment.

## Running

Each compiled binary takes a fixed positional argument list. Run it with **zero** arguments to fall back to
an interactive `cin` prompt that asks for each parameter in turn.

| Driver | Args (after program name) |
|---|---|
| `order_1stocDP_unity.cpp` | `dt t_max L R a_start a_end div dP PREFIX` (9) |
| `order_2stocDP_unity.cpp` | `dt t_max L R a_start a_end div dP init_frac_graz aij_scale PREFIX` (11) |
| `order_2stocDP_burnin.cpp` | ...as above, plus `input_frame_subdir` (12) |
| `order_3stocDP_unity.cpp` | `dt t_max L R a_start a_end div dP init_frac_pred aij_scale ajm_scale PREFIX` (12) |
| `order_3stocDP_burnin.cpp` | ...as above, plus `input_frame_subdir` (13) |

Parameter meanings (same convention as `Rietkerk_FD`):

- `dt` — integration timestep.
- `t_max` — maximum simulated duration.
- `L`/`g` — grid side length.
- `R`/`r` — number of stochastic replicates.
- `a_start`, `a_end` — start/end of the control-parameter (`a`) scan range.
- `div` — number of divisions of the `[a_start, a_end]` range; also doubles as the OpenMP thread count when
  launched from the scheduling scripts.
- `dP` — perturbation/"kick" magnitude used to nudge the system into the high-density state.
- `init_frac_pred` / `init_frac_graz` — fractional deviation of the initial grazer/predator density from
  their mean-field reference value.
- `aij_scale`, `ajm_scale` — scaling factors applied to the base grazer- and predator-attack rates.
- `PREFIX` — string prefix used to name/tag all output files for the run.
- `input_frame_subdir` (burn-in drivers only) — subdirectory to read the pre-burned-in initial frame from.

Example (3-species run, using the binary built above):

```bash
./dp_3sp.out 0.1 400000 128 5 0.04 0.06 8 10000 1.1 1 1 MyRunPrefix
```

Output is written under `../Data/DP/{Frames,Prelims,Stochastic}/<SpB>Sp/<PREFIX>_...`, relative to this
directory — see `/CLAUDE.md` for the full data pipeline (`Reorganised_Frames`, `Data_Processing/frame_builder.py`, etc.).

## Running many jobs (Scheduling Bashes)

For parameter sweeps across many `(a, dt, L, R, ...)` combinations at once, `Scheduling Bashes/` provides
launcher scripts driven by a shared tab-delimited init-file format (sample files under
`Scheduling Bashes/Init_Files/`):

```
dt   t_max   L   R   a_st   a_end   div   dP   PREFIXES
```
One line = one job.

- **`twilight_screening.bash` ** — launch a batch of jobs locally, each in its
  own detached GNU `screen` session, compiling and running the appropriate `order_${SPB}stocDP_*`
  combination per line of the init file.
  ```bash
  ./"Scheduling Bashes/twilight_screening.bash" "Scheduling Bashes/Init_Files/chonKY_init-3SP-DsB6L9-UNITY.txt" 3 1
  # <init_file.txt> [SpB=3] [Init_Type=0]
  ```
  `twilight_screening.bash` links `-larmadillo -lfftw3_threads -lfftw3 -lm`; `midnight_screening.bash` — ensure you have these libraries installed before use. Attach to a running job with `screen -r <name>`; list sessions with `screen -list`.
- **`amarel_job_array_prep.bash`** — generates and submits (`sbatch`) a SLURM job-array script from the same
  init-file format, for running the sweep on the Amarel HPC cluster.
  ```bash
  ./"Scheduling Bashes/amarel_job_array_prep.bash" "Scheduling Bashes/Init_Files/amarel_init-2SP-B6.txt" my_job_name 2
  # <init_file.txt> <job_name> <SPB> [Init_Type=0] [cpus_per_task] [num_jobs]
  ```
  Check the generated `<job_name>_array.sh` before submitting — the FFTW3 library path, SLURM partition,
  memory/time limits, and notification email are hardcoded in the script template (the email is a
  placeholder, `xXx@uni.gov`) and need updating for your account/cluster.
- **`midnight_screening_DDM.bash`** — do not use: it's a stale leftover copy from `Rietkerk_FD` referencing
  files (`rietkerk_bjork_basic-DDM.cpp`, `order_3DDMstoc_test_rietkerk.cpp`) that don't exist in this directory.
