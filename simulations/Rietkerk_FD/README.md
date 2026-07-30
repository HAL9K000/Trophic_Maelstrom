# Rietkerk_FD Documentation

Finite-difference (Dornic-scheme) stochastic implementation of the Rietkerk vegetation–water dryland
model, extended to a 1/2/3-species trophic chain (vegetation → grazer → predator). Simulations run on a 2D
periodic lattice and are integrated in time with a stochastic Runge-Kutta/Dornic scheme, parallelised over
OpenMP threads (one thread per replicate/rainfall-value combination). Species interaction terms ("gamma",
the local grazing/predation coupling) are computed on each timestep either directly over a neighbourhood
(O(N²)) or via a cached FFT convolution (O(N² log N)) using FFTW3, and initial frames can be generated from
closed-form mean-field-theory expressions parsed at runtime with ExprTk. Each run performs a scan over a
range of the rainfall control parameter `a`, and writes CSV frames, per-timestep aggregate ("Prelim") data,
and summary order-parameter statistics to disk.

This document covers the two core implementation files — `rietkerk_bjork_basic.cpp`/`.h` and the
`order_*_unity_rietkerk.cpp`/`order_*_burnin_rietkerk.cpp` drivers that link against them — and how to
compile and run them. The directory also contains a number of variant/scratch files (a locust
wind-dispersal extension, a density-dependent-mortality variant, older test drivers, deprecated
implementations, a CUDA kernel) that are not covered here; consult those files directly if you need them.

For the repository-wide picture (data-output conventions, how this fits with `Percolation_FD`,
`Rietkerk_FPE`, and the `Utilities`/`Data_Processing` post-processing scripts), see the top-level
`simulations` README.

## Dependencies

- `g++` with C++23 support (the scheduling scripts here use `g++-14`)
- OpenMP (`-fopenmp`)
- [FFTW3](http://www.fftw.org/) — FFT-accelerated interaction-kernel convolution
  (`#include <fftw3.h>` in `rietkerk_bjork_basic.h`). Not vendored — build/install it yourself and make sure
  your compiler/linker can find its headers and `libfftw3`/`libfftw3_threads`.
- [ExprTk](https://github.com/ArashPartow/exprtk) (header-only) — symbolic-expression parsing used for
  frame initialization (`#include <exprtk.hpp>`). Not vendored — put it on your include path.
- GNU `screen`, only if using the local parallel job launchers in `Scheduling Bashes/`.
- A SLURM environment (e.g. Rutgers' Amarel cluster), only if using the cluster job-array script.

## Core files

- `rietkerk_bjork_basic.h` / `rietkerk_bjork_basic.cpp` — the shared model implementation, compiled into
  every binary. Provides: grid types over plain nested `std::vector` (`D2Vec_Double`, `D3Vec_Double`,
  `D3Vec_Int`, ...); frame-initialization routines (homogeneous MFT, random speckle, gradient,
  burn-in-from-file, gaussian tear, etc.); the `FFTW3_CentralPlanner` struct that caches FFT plans/kernels
  for the interaction-neighbourhood convolution; the RK4/Dornic stochastic integrators
  (`RK4_Integrate_Stochastic*`, `rietkerk_Dornic_2D*`) for the 1/2/3-species variants; both direct and
  FFT-accelerated gamma-calculation routines (`calc_gamma_*`); and CSV/frame I/O.
- `rietkerk_bjork_constants_{1,2,3}Sp.h` — pulled in by `rietkerk_bjork_basic.h` based on the `-DSPB` value
  at compile time (see below); defines species-count-specific constants (total species count `Sp`, CSV
  frame/prelim column headers).
- `Debug.h` — small third-party (Nadeau Software, CC-BY 3.0) helper for querying peak/current memory (RSS)
  usage, included by `rietkerk_bjork_basic.cpp` for diagnostics.
- `order_{1,2,3}stoc_unity_rietkerk.cpp`, `order_{2,3}stoc_burnin_rietkerk.cpp` — the driver files. Each
  contains only a `main()` and is compiled together with `rietkerk_bjork_basic.cpp` to produce one
  executable. `unity` drivers start from a fresh (random-speckle) initial condition; `burnin` drivers read
  the initial frame from a previously-generated CSV instead.

## Compiling

There is no Makefile/CMakeLists — every run configuration is a standalone `g++` invocation, compiling
`rietkerk_bjork_basic.cpp` together with **one** `order_*` driver file. Which species-count model and which
initial-condition routine gets built is controlled by two **compile-time macros**, not runtime flags. This
mirrors the options exposed by `Scheduling Bashes/twilight_screening.bash`:

```
Usage: twilight_screening.bash <Path/to/init_file.txt> [Optional: <SpB> <Init_Types>]
  SpB:        number of biotic species (default 3)
  Init_Types: 1 = Random MFT-Based Speckles, 2 = Burn-in Frames read from file, 0 = Homogeneous MFT Frames
```

| Macro | Values | Meaning |
|---|---|---|
| `-DSPB=<N>` | `1`, `2`, `3` | Number of biotic species: vegetation only, +grazer, +grazer+predator. Selects which `order_*` driver file to compile against, and which `rietkerk_bjork_constants_<N>Sp.h` gets pulled in. |
| `-DINIT=<N>` | `0`, `1`, `2` | Initial condition: `0` = homogeneous mean-field-theory frame, `1` = random MFT-based speckles ("unity"), `2` = burn-in frame read from a file on disk. `0`/`1` require an `order_*stoc_unity_rietkerk.cpp` driver; `2` requires the corresponding `order_*stoc_burnin_rietkerk.cpp` driver. |

### 1-species (vegetation only)

```bash
g++-14 -O3 -march=native -DSPB=1 -DINIT=1 \
  rietkerk_bjork_basic.cpp order_1stoc_unity_rietkerk.cpp \
  -L${HOME}/.local/lib -lfftw3_threads -lfftw3 -lm -fopenmp \
  -o rietkerk_1sp.out -std=c++23
```

### 2-species (vegetation + grazer)

```bash
g++-14 -O3 -march=native -DSPB=2 -DINIT=1 \
  rietkerk_bjork_basic.cpp order_2stoc_unity_rietkerk.cpp \
  -L${HOME}/.local/lib -lfftw3_threads -lfftw3 -lm -fopenmp \
  -o rietkerk_2sp.out -std=c++23
```

### 3-species (vegetation + grazer + predator)

```bash
g++-14 -O3 -march=native -DSPB=3 -DINIT=1 \
  rietkerk_bjork_basic.cpp order_3stoc_unity_rietkerk.cpp \
  -L${HOME}/.local/lib -lfftw3_threads -lfftw3 -lm -fopenmp \
  -o rietkerk_3sp.out -std=c++23
```

### Burn-in initial condition (`-DINIT=2`)

Swap the `unity` driver for the matching `burnin` one (2- and 3-species only — there is no
`order_1stoc_burnin_rietkerk.cpp`):

```bash
g++-14 -O3 -march=native -DSPB=3 -DINIT=2 \
  rietkerk_bjork_basic.cpp order_3stoc_burnin_rietkerk.cpp \
  -L${HOME}/.local/lib -lfftw3_threads -lfftw3 -lm -fopenmp \
  -o rietkerk_3sp_burnin.out -std=c++23
```

The committed `.vscode/c_cpp_properties.json` is IntelliSense-only and does not reflect these real build
commands — don't treat it as authoritative for compiler/flag choices.

## Running

Each compiled binary takes a fixed positional argument list matching its `argc` check. Run it with **zero**
arguments to fall back to an interactive `cin` prompt that asks for each parameter in turn — a convenient
way to see what each one means without reading the source.

| Driver | Args (after program name) |
|---|---|
| `order_1stoc_unity_rietkerk.cpp` | `dt t_max L R a_start a_end div dP PREFIX` (9) |
| `order_2stoc_unity_rietkerk.cpp` | `dt t_max L R a_start a_end div dP init_frac_graz aij_scale PREFIX` (11) |
| `order_2stoc_burnin_rietkerk.cpp` | ...as above, plus `input_frame_subdir` (12) |
| `order_3stoc_unity_rietkerk.cpp` | `dt t_max L R a_start a_end div dP init_frac_pred aij_scale ajm_scale PREFIX` (12) |
| `order_3stoc_burnin_rietkerk.cpp` | ...as above, plus `input_frame_subdir` (13) |

Parameter meanings:

- `dt` — integration timestep (a value around 0.1 is typical).
- `t_max` — maximum simulated duration (hours).
- `L`/`g` — grid side length (each cell corresponds to `dx = 100 m`).
- `R`/`r` — number of stochastic replicates.
- `a_start`, `a_end` — start/end of the rainfall-parameter (`a`) scan range.
- `div` — number of divisions of the `[a_start, a_end]` range; also doubles as the OpenMP thread count when
  launched from the scheduling scripts.
- `dP` — perturbation/"kick" magnitude used to nudge the system into the high-density state.
- `init_frac_pred` / `init_frac_graz` — fractional deviation of the initial grazer/predator density from
  their mean-field-theory reference value.
- `aij_scale`, `ajm_scale` — scaling factors applied to the base grazer- and predator-attack rates.
- `PREFIX` — string prefix used to name/tag all output files for the run.
- `input_frame_subdir` (burn-in drivers only) — subdirectory to read the pre-burned-in initial frame from.

Example (3-species run, using the binary built above):

```bash
./rietkerk_3sp.out 0.1 400000 128 5 0.04 0.06 8 10000 1.1 1 1 MyRunPrefix
```

Output (CSV frames, "Prelims" aggregate stats, "Stochastic" summary/order-parameter stats) is written under
`../Data/Rietkerk/{Frames,Prelims,Stochastic}/<SpB>Sp/<PREFIX>_...`, relative to this directory — see
`/CLAUDE.md` for the full data pipeline (`Reorganised_Frames`, `Data_Processing/frame_builder.py`, etc.).

## Running many jobs (Scheduling Bashes)

For parameter sweeps across many `(a, dt, L, R, ...)` combinations at once, `Scheduling Bashes/` provides
launcher scripts driven by a shared tab-delimited init-file format (sample files under
`Scheduling Bashes/Init_Files/`):

```
dt   t_max   L   R   a_st   a_end   div   dP   PREFIXES
```
One line = one job.

- **`twilight_screening.bash` ** — launch a batch of jobs locally, each in its own detached GNU `screen` session, compiling and running the appropriate `order_${SPB}stoc_*` combination per line of the init file.
  ```bash
  ./"Scheduling Bashes/twilight_screening.bash" "Scheduling Bashes/Init_Files/chonky_init-REF.txt" 3 1
  # <init_file.txt> [SpB=3] [Init_Type=1]
  ```
  Attach to a running job with `screen -r <name>`; list sessions with `screen -list`.
- **`amarel_job_array_prep.bash`** — generates and submits (`sbatch`) a SLURM job-array script from the same
  init-file format, for running the sweep on the Amarel HPC cluster.
  ```bash
  ./"Scheduling Bashes/amarel_job_array_prep.bash" "Scheduling Bashes/Init_Files/amarel_init-3SP-NREF.txt" my_job_name 3
  # <init_file.txt> <job_name> <SPB> [Init_Type=1] [cpus_per_task] [num_jobs]
  ```
  Check the generated `<job_name>_array.sh` before submitting — the SLURM partition, memory/time limits, and
  notification email are hardcoded in the script template and may need adjusting for your account/cluster.
- Variant scheduling scripts ( `locust_*_screening.bash`,
  `locust_amarel_job_array_prep.bash`, `nonstoc_amarel_job_array_prep.bash`) follow the same pattern for locust-dispersal, and non-stochastic variants respectively — see those files directly.
