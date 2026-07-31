# Rietkerk_FPE Documentation

Fokker-Planck-equation (FPE) variant of the Rietkerk vegetation-water dryland model (see `../Rietkerk_FD` for the base finite-difference model this extends). Instead of treating dispersal as plain diffusion, this variant adds an advection-diffusion step for the population density fields, computed either numerically on the GPU via a dedicated CUDA kernel, or analytically on the CPU via precomputed Gaussian stencils, to capture directional/velocity-biased movement. Simulations run on a 2D periodic lattice, parallelised over OpenMP threads, with the same consumer movement switching behaviour ("gamma") and stochastic-integration machinery as `Rietkerk_FD`. Each run scans a range of the rainfall control parameter `a` and writes CSV frames, per-timestep aggregate ("Prelim") data, and summary order-parameter statistics to disk.

> **Testing status: ⚠️⚠️ ** Only the **CUDA path** (`kernel_advdiff_FKE2D.cu`, `-DBARRACUDA` builds) has been
> rigorously tested and validated. The CPU-only path (compiling without `-DBARRACUDA`, which falls back to
> the analytic Gaussian-stencil integration in `rietkerk_bjork_basic.h`) is more experimental, has seen far
> less validation, and may contain bugs. Use the CUDA build for anything beyond quick local exploration, and
> treat CPU-only FPE results with corresponding caution until cross-checked against the CUDA path.

This document covers the core implementation files — `rietkerk_bjork_FKE.cpp`/`rietkerk_bjork_basic.h`, the
`kernel_advdiff_FKE2D.cu` CUDA kernel, and the `order_*_unity_rietkerk.cpp`/`order_*_burnin_rietkerk.cpp`
drivers that link against them — and how to compile and run them. The directory also contains a
density-dependent-mortality variant, older test drivers, and sample burn-in frames, which are not covered
here.

For the repository-wide picture (data-output conventions, how this fits with `Rietkerk_FD`, `Percolation_FD`,
and the `Utilities`/`Data_Processing` post-processing scripts), see the top-level `simulations` README`.

## Dependencies

- `g++` with C++20/23 support (scheduling scripts here use `g++-12`, `g++-14` in some launchers)
- OpenMP (`-fopenmp`)
- [FFTW3](http://www.fftw.org/) and [ExprTk](https://github.com/ArashPartow/exprtk) (header-only) — same as
  `Rietkerk_FD`; required by the shared `rietkerk_bjork_basic.h` header. Neither is vendored.
- **CUDA toolkit (`nvcc`) — required for the validated path.** `kernel_advdiff_FKE2D.cu` implements the GPU
  advection-diffusion kernel; the CUDA scheduling scripts check for `nvcc` up front and refuse to run
  without it.
- GNU `screen`, only if using the local parallel job launchers in `Scheduling Bashes/`.
- A SLURM environment (e.g. Rutgers' Amarel cluster), only if using the cluster job-array script (see the
  caveat in [Running many jobs](#running-many-jobs-scheduling-bashes) — the one provided here is stale).

## Core files

- `rietkerk_bjork_basic.h` — the same shared header family as `Rietkerk_FD`'s (grid types, frame
  initialization, `FFTW3_CentralPlanner`, RK4/Dornic integrators, gamma calculation), extended with
  FPE-specific machinery: the `MV_INVARIANCE` movement-invariance flag, precomputed Gaussian stencils for
  the analytic distance/velocity FKE, the per-species-pair FKE integrators (`f_DorFKE_2Sp`, `f_DorFKE_3Sp`,
  `advdiff_FKE_MultiSp`), and CUDA kernel-launch declarations gated by `__CUDACC__`/`BARRACUDA`.
- `rietkerk_bjork_constants_{1,2,3}Sp.h` — pulled in by `rietkerk_bjork_basic.h` based on the `-DSPB` value
  at compile time; defines species-count-specific constants.
- `rietkerk_bjork_FKE.cpp` — the FPE-specific implementation (equivalent role to `rietkerk_bjork_basic.cpp`
  in `Rietkerk_FD`), compiled into every binary regardless of whether the CUDA or CPU-only path is used.
- `kernel_advdiff_FKE2D.cu` — CUDA kernel implementing the numerical advection-diffusion solve of the FPE on
  the GPU. This is the code path that has actually been validated; only linked in for `-DBARRACUDA` builds.
- `Debug.h` — small third-party (Nadeau Software, CC-BY 3.0) helper for querying peak/current memory (RSS)
  usage, included by `rietkerk_bjork_FKE.cpp` for diagnostics.
- `order_{2,3}stoc_unity_rietkerk.cpp`, `order_{2,3}stoc_burnin_rietkerk.cpp` — the driver files. Each
  contains only a `main()` and is compiled together with `rietkerk_bjork_FKE.cpp` (and, for CUDA builds,
  linked against `kernel_advdiff_FKE2D.cu`'s object file) to produce one executable. `unity` drivers start
  from a fresh (random-speckle) initial condition; `burnin` drivers read the initial frame from a
  previously-generated CSV instead. There is no 1-species unity/burnin driver in this directory (only a
  scratch `order_1stoc_test_rietkerk.cpp`) — the FPE model here is only exercised for 2- and 3-species runs.

  **NOTE** one FPE-specific behavior: in these drivers, the `dt` command-line argument is accepted but then
  **overridden** with an internally computed `dt = dt_analytical` — the effective integration timestep for
  the FKE step is derived analytically rather than taken from user input.

## Compiling

There is no Makefile/CMakeLists — every run configuration is a standalone `g++`/`nvcc` invocation, compiling
`rietkerk_bjork_FKE.cpp` together with **one** `order_*` driver file (and, for the CUDA path, the compiled
`kernel_advdiff_FKE2D.cu` kernel object). Which species-count model, initial-condition routine, and
integration path get built is controlled by **compile-time macros**, not runtime flags. This mirrors the
options exposed by `Scheduling Bashes/twilight_screening_FPE-CUDA.bash`:

```
Usage: twilight_screening_FPE-CUDA.bash <Path/to/init_file.txt> [Optional: <SpB> <Init_Types> <Mv_Invariance> <ARCH-SM>]
  SpB:           number of biotic species (default 3)
  Init_Types:    type of frame initialisation (default 1)
  Mv_Invariance: 0 = Distance Invariant, 1 = Velocity Invariant (default 0)
  ARCH-SM:       CUDA virtual architecture / streaming multiprocessor (default: left to nvcc)
```

| Macro | Values | Default if omitted | Meaning |
|---|---|---|---|
| `-DSPB=<N>` | `2`, `3` | `3` (`rietkerk_bjork_basic.h` falls back to `#define SPB 3` if omitted) | Number of biotic species. Selects which `order_*` driver file to compile against, and which `rietkerk_bjork_constants_<N>Sp.h` gets pulled in. In practice, always pass this explicitly, matching the driver file you're compiling against. |
| `-DINIT=<N>` | `0`, `1`, `2` | `1` (`#define INIT 1` fallback) | Initial condition: `0` = homogeneous, `1` = random speckles ("unity"), `2` = burn-in frame read from a file on disk. `0`/`1` require an `order_*stoc_unity_rietkerk.cpp` driver; `2` requires the corresponding `order_*stoc_burnin_rietkerk.cpp` driver. |
| `-DMV_INVARIANCE=<N>` | `0`, `1` | `0` (`#define MV_INVARIANCE 0` fallback) | Movement-invariance type for the FKE step: `0` = distance-invariant, `1` = velocity-invariant. |
| `-DBARRACUDA` | (flag) | off | **Enables the validated CUDA/GPU advection-diffusion path.** Requires linking against a compiled `kernel_advdiff_FKE2D.cu` object and the CUDA runtime. Omitting it falls back to the less-tested CPU-only analytic path. |
| `-DDEBUG` | (flag) | off | Enables verbose per-timestep diagnostic `cout`/`cerr` prints throughout `rietkerk_bjork_FKE.cpp` (status updates before/after the advection-diffusion step, `Rho_dt` and gamma sanity checks), plus a safety check on every gamma value: if a computed `gamma` value is NaN or Inf, the current state is saved to an `ERROR_GAMMA_*.csv` frame before the process exits. Without `-DDEBUG`, a NaN/Inf gamma value is only logged, not caught or saved. |

### CUDA build (recommended — the tested path)

Three-stage build: compile the driver and `rietkerk_bjork_FKE.cpp` as CUDA-aware object files, compile the
CUDA kernel with `nvcc`, then link everything together.

2-species:
```bash
g++-14 -DSPB=2 -DINIT=1 -DMV_INVARIANCE=0 -DBARRACUDA -c order_2stoc_unity_rietkerk.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include
g++-14 -DSPB=2 -DINIT=1 -DMV_INVARIANCE=0 -DBARRACUDA -c rietkerk_bjork_FKE.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include
nvcc -O3 --use_fast_math -arch=sm_<XX> -c -o kernel_FKE.o kernel_advdiff_FKE2D.cu -std=c++20
g++-14 -o rietkerk_2sp_cuFKE.out rietkerk_bjork_FKE.o order_2stoc_unity_rietkerk.o kernel_FKE.o \
  -fopenmp -O3 -march=native -std=c++20 \
  -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcudadevrt -Wl,-rpath=/usr/local/cuda/lib64
```

3-species:
```bash
g++-14 -DSPB=3 -DINIT=1 -DMV_INVARIANCE=0 -DBARRACUDA -c order_3stoc_unity_rietkerk.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include
g++-14 -DSPB=3 -DINIT=1 -DMV_INVARIANCE=0 -DBARRACUDA -c rietkerk_bjork_FKE.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include
nvcc -O3 --use_fast_math -arch=sm_<XX> -c -o kernel_FKE.o kernel_advdiff_FKE2D.cu -std=c++20
g++-14 -o rietkerk_3sp_cuFKE.out rietkerk_bjork_FKE.o order_3stoc_unity_rietkerk.o kernel_FKE.o \
  -fopenmp -O3 -march=native -std=c++20 \
  -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcudadevrt -Wl,-rpath=/usr/local/cuda/lib64
```
(replace `sm_<XX>` with your GPU's compute capability, e.g. `sm_89`; omit `-arch=sm_<XX>` to let `nvcc`
choose a default). For the burn-in initial condition (`-DINIT=2`), swap in the matching `*_burnin_*` driver
in place of `*_unity_*` throughout.

### CPU-only build (experimental — use with caution)

```bash
g++-12 -DSPB=3 -DINIT=1 -DMV_INVARIANCE=0 \
  rietkerk_bjork_FKE.cpp order_3stoc_unity_rietkerk.cpp \
  -fopenmp -o rietkerk_3sp_fpe_cpu.out -std=c++23
```
This omits `-DBARRACUDA`, so the advection-diffusion step runs through the analytic Gaussian-stencil code
path in `rietkerk_bjork_basic.h` instead of the GPU kernel — see the testing-status note above.

## Running

Each compiled binary takes a fixed positional argument list. Run it with **zero** arguments to fall back to
an interactive `cin` prompt that asks for each parameter in turn.

| Driver | Args (after program name) |
|---|---|
| `order_2stoc_unity_rietkerk.cpp` | `dt t_max L R a_start a_end div dP init_frac_graz aij_scale PREFIX` (11) — `dt` is accepted but overridden by `dt_analytical` internally |
| `order_2stoc_burnin_rietkerk.cpp` | ...as above, plus `input_frame_subdir` (12) |
| `order_3stoc_unity_rietkerk.cpp` | `dt t_max L R a_start a_end div dP init_frac_pred aij_scale ajm_scale PREFIX` (12) — `dt` overridden as above |
| `order_3stoc_burnin_rietkerk.cpp` | ...as above, plus `input_frame_subdir` (13) |

Parameter meanings match `Rietkerk_FD` (`dt` timestep — see note above, `t_max` max simulated duration,
`L`/`g` grid side length, `R`/`r` replicate count, `a_start`/`a_end` control-parameter scan range, `div`
number of scan divisions / OpenMP thread count, `dP` perturbation kick magnitude, `init_frac_*` fractional
deviation of initial grazer/predator density from the mean-field reference, `aij_scale`/`ajm_scale`
grazer/predator interaction-rate scale factors, `PREFIX` output-file tag, `input_frame_subdir` burn-in
source subdirectory).

Example (3-species CUDA run, using the binary built above):

```bash
./rietkerk_3sp_cuFKE.out 0.1 400000 128 5 0.04 0.06 8 10000 1.1 1 1 MyRunPrefix
```

Output is written under `../Data/Rietkerk/{Frames,Prelims,Stochastic}/<SpB>Sp/<PREFIX>_...`, relative to
this directory — the same convention as `Rietkerk_FD` (CUDA/FPE runs are not segregated into a separate
`Data/` subtree). See `/CLAUDE.md` for the full data pipeline.

## Running many jobs (Scheduling Bashes)

For parameter sweeps across many `(a, dt, L, R, ...)` combinations at once, `Scheduling Bashes/` provides
launcher scripts driven by a shared tab-delimited init-file format (sample files under
`Scheduling Bashes/Init_Files/`):

```
dt   t_max   L   R   a_st   a_end   div   dP   PREFIXES
```
One line = one job. **Naming here is a bit confusing: `twilight_screening_FPE.bash` is not the CPU-only
counterpart to `twilight_screening_FPE-CUDA.bash` — both require and use `nvcc`.** The only genuinely
CPU-only launcher is `midnight_screening_FPE.bash`.

- **`midnight_screening_FPE.bash`** — CPU-only launcher (experimental path). Compiles and runs
  `rietkerk_bjork_FKE.cpp` + `order_${SPB}stoc_{unity,burnin}_rietkerk.cpp` without `-DBARRACUDA`, one job
  per line of the init file, each detached in a `screen` session.
  ```bash
  ./"Scheduling Bashes/midnight_screening_FPE.bash" "Scheduling Bashes/Init_Files/chonKY_init-3SP-ASCALE.txt" 3 1
  # <init_file.txt> [SpB=3] [Init_Type=1]
  ```
- **`twilight_screening_FPE.bash` / `twilight_screening_FPE-CUDA.bash`** — CUDA launchers (recommended,
  tested path), functionally equivalent to each other. Compile with `-DBARRACUDA`, build/link
  `kernel_advdiff_FKE2D.cu`, and run each job detached in a `screen` session.
  ```bash
  ./"Scheduling Bashes/twilight_screening_FPE-CUDA.bash" "Scheduling Bashes/Init_Files/chonKY_init-3SP-ASCALE-HX-CuFKD-SIMP.txt" 3 1 0 89
  # <init_file.txt> [SpB=3] [Init_Type=1] [Mv_Invariance=0] [ARCH-SM]
  ```
  Attach to a running job with `screen -r <name>`; list sessions with `screen -list`.
- **`amarel_job_array_prep.bash`** — **stale/non-functional as committed**: it's a copy-paste leftover from
  `Rietkerk_FD`'s script and compiles `rietkerk_bjork_basic.cpp` (which doesn't exist in this directory) with
  no CUDA step at all. Don't use it without rewriting its compile step to build `rietkerk_bjork_FKE.cpp` +
  the appropriate `order_*` driver, plus the `nvcc`/`-DBARRACUDA` steps if you want the validated GPU path.
- **`midnight_screening_DDM.bash`** — launcher for the density-dependent-mortality variant (out of scope for
  this document; see the file directly if needed).
