# Data_Processing

Visualisation and derived-statistics layer for the reorganised Rietkerk/DP simulation output produced by
`../Utilities/reorganise_dir.py` and `../Utilities/reorganise_prelims_dir.py`. Where `Utilities/` reorganises
raw C++ output and computes per-directory summary statistics (means/std, FFT power spectra, clustering,
potential wells, spatio-temporal synchrony), this directory turns that reorganised data and those derived
statistics into heatmap images, videos, and summary plots.

The active script here is **`frame_builder.py`**. `deprecated_unsupervised_clustering.py` is, per its own
name, deprecated and out of scope for this document.

## Dependencies

- Python 3.10+ (developed/tested against the `$HOME/dask_py3.11env` virtualenv, which has every dependency
  below already installed)
- `numpy`, `pandas`, `scipy`, `regex`, `powerlaw`
- `matplotlib`, `seaborn`, `adjustText` (plotting)
- `opencv-python` (`cv2`, for building videos from the generated PNG frames)
- `scikit-learn` (`sklearn.cluster.KMeans`, used by the clustering-adjacent analysis functions)
- **`ffmpeg`** (external binary, not a Python package) -- used via `subprocess` to re-encode the
  `cv2`-generated `.mp4` videos to H264 (`ffmpeg -c:v libx264 -crf 20 -preset veryslow ...`) for wider
  playback compatibility. This step is Linux-only and is skipped automatically (with a warning) if `ffmpeg`
  isn't found (`shutil.which("ffmpeg")`) or the OS isn't Linux -- video generation itself still works either
  way, you just keep the original (non-H264) `.mp4`.
- Jupyter/`ipykernel`/`nbformat`, only if running `frame_builder_walkthrough.ipynb` (see below)

Run `python3 frame_builder.py --help` for the full list of overridable command-line flags (`--SPB`, `--indir`,
`--outdir`, `--g`, `--dP`, `--Geq`, `--R_max`, `--prefixes`, `--a_vals`, `--T_vals`, `--TS_vals`,
`--TCorr_vals`, `--a_scaling`) with their current hardcoded defaults noted; any flag left unset falls back to
that hardcoded value.

## How this relates to the `Utilities/` directory structure

`frame_builder.py` expects its `in_dir` to already be in the **reorganised** layout that
`Utilities/reorganise_dir.py`/`reorganise_prelims_dir.py` produce (see `Utilities/README.md`), not the raw C++
output layout:
```
in_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv          (Frames)
in_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/TimeSeries/TSERIES_T_{T}_a_{a_val}_R_{R}.csv    (Prelims)
```
It also reads the `a_vals.txt`/`T_vals.txt`/`TimeSeries_vals.txt`/`maxR.txt` manifest files that
`Utilities/reorganise_dir.py`/`reorganise_prelims_dir.py` write into that same tree, and several analysis
functions here consume the *derived* statistics files those scripts' post-processing stages produce in-place
alongside the reorganised data (`MEAN_STD_{FRAME,GAMMA}_*.txt`, `FFT_POWERspectra.csv`,
`CLUST/{binType}_{nbins}/CLUSTERED_FREQUENCIES.csv`, `Pot_Well/{Pot_Well,LOCAL_MINIMA}.csv`,
`2DCorr/{Auto,Cross}_*.csv`, `2DCorr/FFT/HarmonicPeaks_*.csv`). If a given post-processing stage was never run
over a dataset, the corresponding `frame_builder.py` analysis simply has nothing to read (see
[Known edge cases and bugs](#known-edge-cases-and-bugs) below for how that's handled -- or isn't).

`frame_builder.py` writes everything it generates -- PNG frames, videos, summary plots, and the `DEBUG_*.csv`
tables some functions produce -- to a separate `out_dir`, mirroring `in_dir`'s directory structure underneath
it, rather than writing back into `in_dir`.

## Functionality overview

`frame_builder.py` is organised around two closely related data sources, matching the two `Utilities/`
reorganisation scripts:

- **Frame data** (per-cell spatial snapshots, one row per grid cell, columns `a_c, x, <species...>`):
  - `frame_visualiser` renders each snapshot as a heatmap PNG (per species, plus a combined subplot); `home_video`
    stitches a directory of such PNGs into an `.mp4`. Together these are what most of the driver section's video
    blocks use, in several groupings (all replicates combined, one video per replicate, across-`a` at fixed `T`,
    with a GAMMA overlay for spreading-test runs).
  - `analyse_FRAME_EQdata` plots the equilibrium (long-time-averaged) state vs. the control parameter `a`.
  - `analyse_FRAME_FFTPOWERdata` plots each replicate's 2D spatial FFT power spectrum as a heatmap (frequency
    vs. `a`), for spotting a characteristic pattern wavelength.
  - `analyse_FRAME_POTdata` plots KDE-based potential wells (and optionally their local minima) vs. `a`, for
    visualising bistability/alternative-stable-states behavior; also usable on the GAMMA/interaction fields.
  - `analyse_FRAME_CLUSTEREDISTdata` and `analyse_FRAME_CLUSTERSTATS_timeseriesData` plot cluster-size
    frequency distributions and macro cluster statistics (mean cluster size, count, density, ...) vs. `a`.
  - `analyse_FRAME_CORRdata` and `analyse_FRAME_FFTCORRdata` plot the spatio-temporal
    synchrony/cross-correlation data (NCC/ZNCC/AMI/MI/bivariate Moran's I) and the harmonic frequencies
    extracted from those correlation curves.
- **Prelims data** (spatially-pooled time series, one row per timepoint):
  - `analyse_PRELIMS_TIMESERIESdata` plots each species' time series, with optional power-law/decay fit
    overlays; also used (with different labels) for the movement/GAMMA (`"MOVSERIES"`) Prelims files.
  - `analyse_PRELIMS_EQdata` plots the equilibrium (time-window-averaged) state vs. `a`, optionally as violin
    plots of the per-replicate distribution, and writes the `DEBUG_*.csv` tables that...
  - `multiplot_PRELIMS_CompareEQ` later combines across multiple prefixes / grid sizes / `dP` / `Geq` values
    into a single comparison plot.
  - `analyse_PRELIMS_TRAJECTORYdata` plots phase-space trajectories (e.g. grazer vs. predator concentration
    over time), per replicate and combined.

Species count (`SPB`) determines which state-variable/movement-label lists (`variable_labels`,
`var_frame_labels`, `move_var_labels`) and default time-averaging windows get used throughout -- see the
walkthrough notebook for the exact `SPB == {1,2,3}` derivation.

**For a full, runnable walkthrough of every one of the above functions** -- including the several
call-pattern variants that exist only as commented-out alternatives in `frame_builder.py`'s own driver
section, each with a markdown description of what it does and what it expects on disk -- see
**`frame_builder_walkthrough.ipynb`** in this directory. It's built to be run cell-by-cell against your own
`in_dir`/`out_dir`/`prefixes`/etc. (redefine them at the top of the notebook; see the notebook's own note
about the handful of globals -- `g`, `dP`, `Geq`, `SPB`, `R_max`, `in_dir`, `out_dir` -- that must be set via
`fb.<name> = ...` rather than a plain local variable, since some functions read them as bare module globals
rather than as parameters).

## Known edge cases and bugs

- `analyse_PRELIMS_TIMESERIESdata` blocks on an internal `input("Press F to Continue...")` prompt partway
  through its run (after it has already generated its plots and video) -- not a bug, but worth knowing before
  running it non-interactively (e.g. from a plain script or `nbconvert --execute`), where it will raise
  `EOFError` there instead of proceeding. In a normal Jupyter cell this just appears as an input box.
