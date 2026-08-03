# Utilities

A CPU/GPU-agnostic Python (plus a Bash/`rsync`/`expect` component) data pipeline for collating and
post-processing raw simulation output from the C++ models (`../Rietkerk_FD`, `../Rietkerk_FPE`, `../Percolation_FD`). Simulations are frequently run across several remote machines/HPC clusters, so this layer's job starts before reorganisation even begins: `dreamcatcher_automata.bash` can `rsync` each remote device's raw output back to a single host (see [Batch orchestration](#batch-orchestration-dreamcatcher_automatabash)
below), which then get reorganised by the same pipeline into a predictable, analysis-ready structure, with a battery of spatial and temporal statistics computed over them along the way — replicate-averaging, spatial pattern (cluster/FFT/potential-well) analysis, and spatio-temporal synchrony analysis between snapshots.

The two scripts that matter most are **`reorganise_dir.py`** (per-cell "Frame" snapshots) and
**`reorganise_prelims_dir.py`** (aggregate "Prelims" time series) — this document is built around them. Both import their actual numerical machinery from **`glow_up.py`** (`from glow_up import *`), which in turn routes its array operations through **`slick.py`**, [a thin CPU/GPU abstraction layer](https://github.com/HAL9K000/slick). This document also covers **`synthetic_clump_generator.py`** (synthetic initial-condition frame generation), **`gmm_classifier.py`**
(the Python side of the live C++↔Python GMM-clustering bridge used by `Percolation_FD`), and the batch
orchestration script **`dreamcatcher_automata.bash`** (which in turn uses **`expect_commands.sh`** for the
`rsync`-based remote-device collation step).

Everything else in this directory has been consolidated out of the way and is out of scope for this
document: `LegacyUtils/` (`legacy_glow_up.py`, `legacy_reorganise_dir.py`, `legacy_reorganise_prelims_dir.py`,
`CPU_fallback_legacy.py`), `Deprecated_Assets/` (`deprecated_functions.py`, and the now-deprecated
`prelims_automata.bash` — see [Batch orchestration](#batch-orchestration-dreamcatcher_automatabash) below), `DirRenamingUtils/` (`scient_literal_rename.py`, `simple_dir_rename.py`), plus `zipper_merge.bash`.

For the repository-wide picture (how this fits with the C++ simulators and `Data_Processing/`), see the top-level `simulations` README.

## Dependencies

**Setting up dependencies.** `../environment.yml` and `../pip3-requirements.txt` cover everything
listed below (both core and optional) for this layer *and* `Data_Processing/`, since the two share
one environment in practice. Pick whichever of pip or conda you prefer:

- **pip** (venv named `trophic-maelstrom`, matching the convention used elsewhere in this repo):
  ```bash
  python3 -m venv ~/trophic-maelstrom
  source ~/trophic-maelstrom/bin/activate
  pip3 install --upgrade pip
  pip3 install -r ../../pip3-requirements.txt
  ```
  or automate the above, including CUDA-version detection for GPU support, with
  `../setup_environment.sh`.
- **conda/mamba** (environment named `trophic-maelstrom`):
  ```bash
  conda env create -n trophic-maelstrom -f ../../environment.yml
  conda activate trophic-maelstrom
  ```

- Python 3.10+ (Note that the automation scripts invoke `python3.11` explicitly)
- `numpy`, `pandas`, `scipy`, `regex`
- `scikit-learn` (`sklearn.cluster.KMeans`, `sklearn.mixture.GaussianMixture`)
- `scikit-image`, `esda`, `libpysal` (spatial statistics — Moran's I)
- `dask`/`dask.distributed` and/or `joblib` (parallel post-processing) — NOTE: Dask is **strongly recommended** over joblib for concurrent CPU multi--threading over shared GPU devices, as it negates Python's GIL.
- Optionally: `cupy` + `cupyx.scipy.{fft,signal,ndimage}` for GPU acceleration (see [GPU vs
  CPU](#gpu-vs-cpu-when-to-use---gpu)) — entirely optional but **recommended** for grid sizes `L² >= 256x256`, everything degrades gracefully to NumPy/SciPy if
  `cupy` isn't importable or `USE_GPU`/`--gpu` isn't set.
- `7z` (`p7zip`), only if using `dreamcatcher_automata.bash`'s compression step.
- `rsync` and `expect`, only if using `dreamcatcher_automata.bash`'s remote-device data-collation step (see
  below) — `expect` in particular is required by `expect_commands.sh`, which drives the `rsync` calls.

Both `reorganise_dir.py --help` and `reorganise_prelims_dir.py --help` will print argparse's own
authoritative, live description of every flag — run that whenever this document and the script's actual
behavior seem to disagree; the flag descriptions below are a curated overview, not a substitute.

## Directory structure: original → reorganised

Both scripts expect the **original** C++ output layout and reorganise it into the same **new** layout,
differing only in whether they're handling per-cell frame snapshots or aggregate time-series files.

### Frames (`reorganise_dir.py`)

Original (as written by the C++ simulators, e.g. under `../Data/Rietkerk/Frames/Stochastic/<SpB>Sp/`):
```
root_dir/{PREFIX}*_dP_{dP}_Geq_{Geq}/FRAME_*_G_{g}_T_{T}_*_a_{a_val}_*_[jID_{jID}_]*_R_{R}.csv
root_dir/{PREFIX}*_dP_{dP}_Veq_{Veq}/FRAME_*_G_{g}_T_{T}_*_a_{a_val}_*_[jID_{jID}_]*_R_{R}.csv
```
Reorganised:
```
out_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv
out_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Veq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv
```
(`GAMMA_...` instead of `FRAME_...` for gamma/interaction-field files — detected by `"GAMMA"` appearing
anywhere in the source filename.)

### Prelims (`reorganise_prelims_dir.py`)

Original (under `../Data/Rietkerk/Prelims/Stochastic/<SpB>Sp/`):
```
root_dir/{PREFIX}*_dP_{dP}_Geq_{Geq}/TimeSeries/PRELIM_TSERIES_*_G_{g}_T_{T}_*_a_{a_val}_*_R_{R}.csv
root_dir/{PREFIX}*_dP_{dP}_Geq_{Geq}/TimeSeries/PRELIM_MOVSERIES_*_G_{g}_T_{T}_*_a_{a_val}_*_R_{R}.csv
```
Reorganised:
```
out_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/TimeSeries/TSERIES_T_{T}_a_{a_val}_R_{R}.csv
out_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/TimeSeries/MOVSERIES_T_{T}_a_{a_val}_R_{R}.csv
```
(the retained filename prefix, `TSERIES` or `MOVSERIES`, is taken directly from the second underscore-token
of the source filename, so it tracks whatever the C++ side actually wrote.)

### Required input filename format

Both scripts **assume** source filenames contain, in this order:

1. `_G_{L}_T_{Time}_` — grid size and the frame's/snapshot's time value,
2. `_a_{a_val}_` — the scan-parameter value,
3. `_R_{R}.csv` — the replicate index,

with `_jID_{jID}_` (a per-replicate job ID assigned by the HPC scheduler) optionally present anywhere in
between. **Filenames that don't follow this pattern will not be picked up** — this is a hard requirement,
not a best-effort heuristic; there is no fallback parser.

### Copy/rename engine and collision handling

Both scripts' `main()` function walks `root_dir` for subdirectories matching the requested
`{PREFIX}*_dP_{dP}_Geq/Veq_{val}` pattern, discovers the unique `a`, `T` (and `jID`, if present) values
present, and copies+renames matching files into the new tree with `shutil.copy2` (preserving metadata). If
the destination file already exists:

- If it's literally the same file (checked by inode, via `os.fstat().st_ino`/`os.path.samefile`), the copy is
  silently skipped.
- Otherwise, it is treated as a **distinct replicate** — the incoming file's `R` index is bumped to
  `max_R + 1` among existing files of that type in the destination directory, so nothing is ever silently
  overwritten.

A `savedir.txt` file is appended to each source subdirectory recording which output directories it has
already been reorganised into, so re-running the script against the same source data doesn't reprocess
directories it already handled (you'll get an interactive/logged warning if it detects this). `T_vals.txt`
(Frames) / `TimeSeries_vals.txt` (Prelims) and `a_vals.txt` manifest files are written into the output tree
to record which values are actually present, and later stages of the pipeline (and `Data_Processing/frame_builder.py`)
read these back rather than re-deriving them from disk listings every time.

## `reorganise_dir.py`

```
python3.XX reorganise_dir.py [--dynamic] [--gpu] [--CPUCores N] [--prefixes P [P ...]]
                              [--indir DIR] [--outdir DIR] [--dP N]
                              [--Geq VAL|NA] [--Veq VAL|NA] [--L N [N ...]]
                              [--indx_vals_t N] [--tmin N] [--tmax N]
```
Run `python3.XX reorganise_dir.py --help` for argparse's own descriptions. Key flags:

| Flag | Meaning |
|---|---|
| `--indir` | Root directory holding the original C++ output. |
| `--outdir` | Root of the reorganised output tree. |
| `--prefixes` | One or more `{PREFIX}` values to process (each becomes a top-level subdirectory of `outdir`). |
| `--dP` | The `dP` (perturbation kick) value embedded in the source subdirectory names. |
| `--Geq` / `--Veq` | The `Geq`/`Veq` value embedded in the source subdirectory names; pass `NA` for whichever one isn't used. |
| `--L` | One or more grid sizes to process. |
| `--indx_vals_t` | How many `T` values to extract per `a`: `-n` = the `n` **largest** (latest) T values, `n` = the `n` **smallest** (earliest). |
| `--tmin` / `--tmax` | Restrict extracted `T` values to this range (in addition to `--indx_vals_t`). |
| `--CPUCores` | Worker count for the parallelised post-processing stage (clamped to `os.cpu_count()`). |
| `--gpu` | Route the numerically-heavy post-processing (especially synchrony analysis) through `cupy` via `slick.py`, if available. |
| `--dynamic` | Skip flag parsing entirely and instead interactively prompt for every parameter at the terminal (useful for one-off exploratory runs). |

Once `main()` has finished reorganising files, the script runs a post-processing pass over the new tree. As
shipped, the `__main__` block runs `main()` then `unified_post_process(..., "FRAME", ...)` and
`unified_post_process(..., "GAMMA", ...)`; the synchrony-analysis call (`post_imgprocess`) and the
second-stage harmonic-frequency pass (`post_process_df_files`) are present in the file but need to be
explicitly enabled — see the `__main__` block if you want them running automatically on every invocation.

### Per-snapshot / per-directory statistics

Run per output subdirectory (i.e. per unique `L`, `a`, `dP`, `Geq/Veq` combination), across all replicates
and all extracted `T` values found there:

- **Replicate mean/std/survival counts** (`gen_MEAN_SD_COLSfiledata`, written as
  `MEAN_STD_{FRAME|GAMMA}_Surviving_Runs.txt` / `..._All_Runs.txt`) — for each state variable column, the
  mean and standard deviation across replicates, computed both over *all* replicates and over only
  *surviving* (non-extinct) replicates, plus the surviving-replicate count at each `T`.
- **Per-cell replicate means** (`gen_MEAN_INDVL_Colsfiledata`, FRAME files only, written as
  `MEAN_FRAME_REPLICATES.txt`) — the replicate-averaged value at *every individual grid cell*, rather than a
  single spatially-pooled statistic; this is what you want for reconstructing a "typical" spatial pattern
  frame rather than a scalar summary.
- **2D spatial FFT power spectra** (`gen_FFT_PowerSpectra`, FRAME files only, written to
  `FFT_PowerSpect/FFT_POWERspectra.csv`) — the 2D Fourier power spectrum of each replicate's spatial field,
  radially binned (`bin_mask="GMM"` by default in `unified_post_process` — bin edges are chosen via a
  Gaussian-mixture fit rather than a fixed bin width, `binwidth` sets the fallback linear bin width). This is
  how you detect a characteristic pattern wavelength (e.g. banding or spot spacing) from a snapshot without
  visually inspecting it.
- **Clustering** (`gen_clustered_data`, written to `CLUST/{cluster_type}_{n_clusters}/{CLUSTERED_FREQUENCIES,MACRO_CLUSTERSTATS,MACRO_FRAMESTATS}.csv`)
  — segments each column's spatial field into `n_clusters` intensity classes and then runs periodic-boundary
  connected-component labelling (`scipy.ndimage.label`) on the resulting binary mask to identify discrete
  spatial clusters/patches. Reports, per column and per replicate: mean cluster size, cluster-size variance,
  number of clusters, cluster density, occupied-site fraction, and mean squared distance from cluster centers
  to both the grid center and the cluster's own center of mass. The classification method is chosen with
  `cluster_type`:
  - **`"KMeans"`** (the default) — fast, deterministic (fixed `random_state=42`), splits values into
    `n_clusters` groups by proximity in 1D intensity space (`sklearn.cluster.KMeans`, `k-means++` init).
    Good default choice for well-separated bimodal (high/low density) fields.
  - **`"GMM"`** — fits a `sklearn.mixture.GaussianMixture` (full covariance, `n_components=n_clusters`)
    instead; softer/probabilistic boundary between classes, more appropriate when the intensity distribution
    isn't cleanly bimodal or has overlapping tails. This is the same clustering approach used by
    `gmm_classifier.py` (see [below](#gmm_classifierpy--cpythonnumpy-c-api-bridge)), just applied here to
    already-written CSV frames instead of live in-memory C++ arrays.
  - **`"Zero"`** — a cheap non-ML fallback: any strictly-positive value is cluster 1, everything else is
    cluster 0. Useful when you already know the field is a clean presence/absence mask and want to skip
    fitting a model entirely.
  Cluster labels are always re-ordered after fitting so that cluster `0` corresponds to the class with the
  lowest center, regardless of which algorithm produced them, so downstream column names stay consistent
  across `cluster_type` choices.
- **Potential-well / KDE analysis** (`gen_potential_well_data`, written to
  `Pot_Well/{Pot_Well,LOCAL_MINIMA}.csv`) — pools each column's values across all replicates at a given `T`
  and estimates a kernel density estimate (KDE) over that pooled distribution, then (optionally) locates its
  local minima. This is the standard way of visualising bistability/alternative-stable-states behavior in
  these models: a KDE with two peaks separated by a well corresponds to a system sitting near a
  degraded/vegetated (or extinct/surviving) bistable transition, and the local minima mark the unstable
  "watershed" density between the two basins.

### Spatio-temporal synchrony/correlation analysis (`post_imgprocess`)

The distinguishing capability of `reorganise_dir.py` over `reorganise_prelims_dir.py`. For every pair of
extracted `T` values within a directory (auto-correlation: a field against itself at a different time;
cross-correlation: different state-variable columns against each other, optionally across different times
too), `post_imgprocess` calls `gen_2DCorr_data` to compute, per pair:

- **NCC / ZNCC** — 2D FFT-based normalized cross-correlation and its zero-mean ("Z") variant (equivalent to a
  spatial Pearson correlation), reporting both the peak correlation value and the spatial lag/index at which
  it occurs.
- **AMI / MI** — adjusted and raw mutual information between the two fields, capturing nonlinear
  dependence that a linear correlation measure like NCC would miss.
- **Bivariate Moran's I** (`calc_Morans`, off by default — enable explicitly if wanted) — a spatial
  autocorrelation statistic (via `esda`/`libpysal`, CPU-only regardless of `--gpu`) measuring whether high/low
  values in one field cluster near high/low values in the other across space.

Results are written per-timepoint-pair to `T_{Tmax}/2DCorr/{Auto,Cross}_{NCC,ZNCC,AMI,MI,BVMoransI}_TD_{delay}.csv`
and pooled across all timepoint pairs in a directory to `2DCorr/{Auto,Cross}_..._T0_{min}_T1_{max}.csv`. A
second pass, `post_process_df_files` (feeding into `get_1D_HarmonicFreq_Prelimsdata`), then runs a 1D FFT over
each of these lag-correlation curves to extract dominant oscillation frequencies/periods and their peaks
(`maxima_finder="find_peaks"` by default; `"cubic_spline"` is also available for smoother interpolated peak
localisation), written to `2DCorr/FFT/{FFTSig,HarmonicPeaks}_..._TD_{delay}.csv` — this is how you'd quantify,
e.g., a periodic predator-prey oscillation's actual period directly from the correlation structure rather than eyeballing a time series.

### GPU vs. CPU: when to use `--gpu`

The synchrony analysis above (2D FFT cross-correlation over every timepoint pair, for every column pair, for every replicate) is by far the most computationally expensive stage in this pipeline. ⚠️⚠️ **For grid sizes `L>= 256`, running this on a GPU is substantially faster than the CPU-only path** — use `--gpu` whenever you're running `post_imgprocess`/synchrony analysis at these grid sizes on a machine with a usable CUDA GPU and `cupy` installed. Also note, combining the `--gpu` path with `dask` as opposed to `joblib` is **strongly recommended for concurrent CPU-GPU interops on shared GPU resources**, as the former avoids the Python GIL for each CPU thread. For smaller grids, or for the copy/rename and basic per-snapshot statistics stages, the CPU path is generally fine and simpler to reason about (no GPU memory/oversubscription concerns — see the [warning below](#gpu-oversubscription-warning) about combining `--gpu` with batch automation). For further details on how GPU--CPU interops is achieved, refer to the [slick documentation](https://github.com/HAL9K000/slick). 

## `reorganise_prelims_dir.py`

Structurally the Prelims analogue of `reorganise_dir.py` — same subdirectory discovery, collision handling,
and manifest-file bookkeeping in `main()` — but with two differences:

```
python3.XX reorganise_prelims_dir.py [--default] [--dynamic] [--CPUCores N] [--prefixes P [P ...]]
                                      [--indir DIR] [--outdir DIR] [--dP N]
                                      [--Geq VAL|NA] [--Veq VAL|NA] [--L N [N ...]]
                                      [--indx_vals_t N] [--dt FLOAT]
```
Run `python3.XX reorganise_prelims_dir.py --help` for the live, authoritative flag list.

1. There is an additional required-in-practice flag, **`--dt`**, since Prelims time series need a known
   sampling interval to be resampled/aligned (`gen_MEAN_INDVL_Prelimsfiledata(..., dt=dt, handle_nonstandardtime_fileconflicts="interpolate")`
   — replicates whose recorded timepoints don't line up exactly are handled by interpolation, not dropped).
   There's also a `--default` flag (distinct from `--dynamic`) that just uses the script's hard-coded default
   values without prompting.
2. **No synchrony/correlation analysis is provided here at all** — there is no equivalent of
   `post_imgprocess`. `post_process()` computes, for each `T`: replicate mean/variance/survival statistics
   over the `TSERIES_*` files (all columns, written to `MEAN_TSERIES_T_{t}.csv`) and over the `MOVSERIES_*`
   files (movement/gamma-related columns only — `<GAM[...]>_x`, `<vx[...]>_x`, `<vy[...]>_x` per species,
   written to `MEAN_MOVSERIES_T_{t}.csv`), plus `maxR_T_{t}.txt`/`movmaxR_T_{t}.txt` replicate-count
   manifests. Both `main()` and `post_process()` run unconditionally from `__main__` (unlike
   `reorganise_dir.py`, nothing here needs manually uncommenting to run by default).

## Batch orchestration (`dreamcatcher_automata.bash`)

For running the same reorganisation+post-processing pipeline over many `(indir, outdir, dP, Geq/Veq, L,
prefixes, ...)` combinations without hand-invoking Python repeatedly, `Automata_Inputs/` holds sample input
files where **each line is a literal argparse flag string**, e.g.:
```
--indir ../Data/Rietkerk/Frames/Stochastic/2Sp/  --outdir ../Data/Rietkerk/Reorganised_Frames/Stoc/2Sp/StdParam_20_MFT/  --dP 10000  --Geq NA  --Veq 7.4774  --indx_vals_t -25  --tmin 60000 --tmax None --L 128 --prefixes DiC-S8LI DiC-S7LI
```

### `dreamcatcher_automata.bash`

```
./dreamcatcher_automata.bash <in_reorganise_dir.txt> <in_reorganise_prelims_dir.txt> <time_gap_hours> <stop_iter> [optionalcommands.txt <ncores>]
```
For `stop_iter` iterations, spaced `time_gap_hours` apart:

1. **(Optional) collate raw output from remote devices.** If an `optionalcommands.txt` file is passed,
   `dreamcatcher_automata.bash` sources `expect_commands.sh` and calls its `execute_commands()` function on
   that file at the start of every iteration. The file is a flat list of `#CMD`/`#INPUT` pairs — one shell
   command per `#CMD` line (in practice, an `rsync` pulling a remote host's `simulations/Data/` tree down into
   a local subdirectory), followed by an `#INPUT` line giving the text to auto-supply if that command prompts
   for input (e.g. an SSH password/passphrase). `execute_commands()` turns each pair into a small
   auto-generated `expect` script (`spawn <command>; expect "*password*" { send "<input>\r" }`) and runs it,
   so a whole batch of remote `rsync` pulls can run unattended. Sample file,
   `Automata_Inputs/addn_commands.txt`:
   ```
   #CMD
   rsync -trlpzv <remotedir1> <localdir1>
   #INPUT
   <placeholder if using key-based auth with none needed>
   #CMD
   rsync -trlpzv --exclude-from='Automata_Inputs/addn_comm_exclude.txt' <remotedir2> <localdir2>
   #INPUT
   <placeholder if using key-based auth with none needed>
   ```
   This is the step that turns "one dataset per remote HPC cluster/device into "one local `Data/` tree", before any of the reorganisation below runs against it. While possible, storing real SSH credentials in a plaintext `#INPUT` line is **strongly discouraged**, and prefer key-based auth with an empty/placeholder `#INPUT` where possible.
2. Reads `in_reorganise_dir.txt` **line by line**; for each line, if the line's `--outdir` already exists
   it's renamed aside with an `HHMM` timestamp suffix (so a repeated run never silently merges into old
   output), then runs `python3.11 reorganise_dir.py <that line> [--CPUCores ncores]` — the entire line,
   **including any `--gpu` flag you put in it**, is passed straight through.
3. Does the same for `in_reorganise_prelims_dir.txt` → `reorganise_prelims_dir.py`.
4. Compresses every output directory touched this iteration with `7z`.
5. Logs errors to `error_log.txt` and sleeps `time_gap_hours` before the next iteration.

Each line is run as its own **sequential, single, foreground `python3.11` process** — there's no
line-level parallelism inside this script itself.

#### GPU oversubscription warning

> **WARNING: ⚠️⚠️ Generally avoid using `--gpu` with `dreamcatcher_automata.bash`** for unattended/batch runs.
> Note `slick.py` supports device-selection policies (round-robin, `least-busy` — picking whichever visible device
> currently has the most free memory, `explicit`), which is correctly dispatched to every Dask worker via
> `daskclient.run(gpu.setup_daskworker_gpu_context, policy='least-busy')` in `unified_post_process` and `gen_2DCorr_data`.
> However, the caution remains because `least-busy` selection is a point-in-time free-memory check, not a lock or
> reservation — concurrent processes (multiple `dreamcatcher_automata.bash` instances run in separate `tmux`
> panes, or a GPU-flagged line combined with a high `--CPUCores`) can each independently see the same device
> as "least busy" and all land on it at once, since nothing in this pipeline coordinates GPU usage *across*
> separate processes. For GPU-accelerated synchrony analysis, prefer invoking `reorganise_dir.py --gpu`
> directly, one job per GPU at a time (and 10-20 CPU threads per GPU), so you retain control over what's actually running concurrently;
> reserve `dreamcatcher_automata.bash` for CPU-only batch runs, or for GPU runs where you've confirmed only
> one `--gpu` invocation will ever be in flight at once.

### `prelims_automata.bash` (deprecated)

The earlier, much thinner single-pass Prelims-only batch runner has been retired to `Deprecated_Assets/` — `dreamcatcher_automata.bash` above is now the one script for both Frames and Prelims batch runs (with
iteration, output-directory collision handling, compression, and remote-device collation that the older script never had).

## `synthetic_clump_generator.py` and the `Input/` file format

Generates synthetic initial-condition frames for seeding new simulation runs, rather than post-processing
existing output. Its main pattern ("HEXBLADE") tiles the grid with a square/simplified hexagonal-center
lattice (true hexagon geometry is present in the code but currently disabled in favour of the simplified
square tiling), places non-periodic 2D Gaussian "clumps" at a configurable random subset of those centers
(`percent_missing` fraction dropped, or a fixed `retain_cluster` count kept), and supports an optional value
gradient across the field (`gradient_type` ∈ `{"linear_X", "linear_-X", "linear_Y", "linear_-Y", "radial",
"Random", None}`). There is no CLI/argparse here — parameters (`SPB`, `L`, replicate count `r`,
`hex_length`, Gaussian radius/amplitude per species, etc.) are set as module-level globals; edit the top of
the file directly for a new configuration.

Output is written to `../Input/Rietkerk/{SPB}Sp/{prefix}/L_{L}_a_0/SEP_{sep}_WID_{wid}/[{grad}_GR/]MSF_{missing}/MAMP_{mamp}_MI(N)VAL_{...}/FRAME_T_0_a_0_R_{replicate}.csv`,
with the same collision-avoidance (bump the replicate index rather than overwrite) as the reorganisation
scripts.

**Input file format**: these generated files, and the `Input/` directory more generally, use exactly the
same per-cell CSV schema as the C++ simulators' own `FRAME_*.csv` snapshot output — a header row of `a_c, x,
<species columns...>` followed by one row per grid cell (`L*L` rows total, `x` the flattened row-major cell
index, `a_c` fixed at `0` for a `T=0` seed frame, and one column per state variable, e.g. `  P(x; t)` for the
plant-biomass field — note the source data really does have leading whitespace baked into that column name,
so preserve it if you're consuming these files programmatically). Because the schema matches the simulators'
own frame format exactly, files under `Input/` can be pointed at directly as burn-in/initial-condition input
(see the `input_frame_subdir` argument documented in the `Rietkerk_FD`/`Percolation_FD` READMEs) without any
conversion step.

## `gmm_classifier.py` — CPython/NumPy C API bridge

Not a standalone script — it's imported live by the C++ simulation code (`Percolation_FD/py_GMM_embedder.hpp`,
active under `-DCRITGMM` builds; see `../Percolation_FD/README.md`) via an embedded CPython interpreter, so
that spatial clustering statistics can be computed **in-memory, during a running simulation**, without a
round-trip through disk. It's deliberately dependency-light — only `numpy`, `sklearn.mixture.GaussianMixture`,
and `scipy.ndimage` — so it can be embedded cheaply.

Its single entry point, `gen_clustered(Rho_t, const_index, L, aval, Tval, Rval, n_clusters=2, periodic=True,
evaluate_clusterfrequencies=True, gen_GMMframe=False)`, takes a raw `(Sp, L*L)` density array straight from
the C++ side's in-memory state and returns scalar cluster statistics (density, occupied-site fraction, mean
cluster size and its variance, number of clusters, and mean squared distance from both the grid center and
each cluster's own center of mass), an optional per-species cluster size/count/frequency table, and an
optional binary cluster-mask frame.

For internal reference: this is algorithmically **the same routine as `gen_clustered_data` in
`glow_up.py`** (specifically its `cluster_type="GMM"` path) — normalize each species field by its column
max, fit a `GaussianMixture`, reorder labels so cluster 0 is the lowest-mean class, binarize, then run
`scipy.ndimage.label` with an explicit periodic-boundary label-stitching pass (merging clusters that wrap
across the grid edges) to get final connected-component clusters — just reimplemented standalone (no
`glow_up`/`slick` import) so it can run inside an embedded interpreter with a minimal dependency footprint, against in-memory C++ arrays instead of CSV files on disk.
