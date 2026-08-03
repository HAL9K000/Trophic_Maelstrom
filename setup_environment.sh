#!/bin/bash
# Sets up a pip venv for the Utilities/ and Data_Processing/ Python layers from
# pip3-requirements.txt, then optionally installs a matching `cupy` build for GPU acceleration
# (skipped automatically if no CUDA toolkit/driver is detected -- everything in both layers works
# fine on CPU alone). If you'd rather use conda/mamba, see environment.yml instead -- this script
# is the pip/venv path only.
#
# Usage:
#   ./setup_environment.sh                          # create ~/trophic-maelstrom, install deps,
#                                                     # auto-detect/install cupy if CUDA is found
#   ./setup_environment.sh --no-gpu                  # skip the cupy step entirely
#   ./setup_environment.sh --venv-path ~/some/where  # use a different venv location

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REQ_FILE="$SCRIPT_DIR/pip3-requirements.txt"
VENV_PATH="$HOME/trophic-maelstrom"
SKIP_GPU=false

while [ $# -gt 0 ]; do
    case "$1" in
        --no-gpu) SKIP_GPU=true; shift ;;
        --venv-path) VENV_PATH="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--no-gpu] [--venv-path PATH]"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [ ! -f "$REQ_FILE" ]; then
    echo "Error: could not find pip3-requirements.txt at $REQ_FILE"; exit 1
fi

if ! command -v python3 &> /dev/null; then
    echo "Error: python3 not found on PATH."; exit 1
fi

if [ -d "$VENV_PATH" ]; then
    echo "Venv already exists at $VENV_PATH -- reusing it and (re)installing requirements into it."
else
    echo "Creating venv at $VENV_PATH ..."
    python3 -m venv "$VENV_PATH"
fi

VENV_PIP="$VENV_PATH/bin/pip3"
VENV_PYTHON="$VENV_PATH/bin/python3"

echo "Upgrading pip in the venv ..."
"$VENV_PIP" install --upgrade pip

echo "Installing dependencies from $REQ_FILE ..."
"$VENV_PIP" install -r "$REQ_FILE"

# --- Optional GPU step: detect CUDA and install a matching cupy build -----------------------
if [ "$SKIP_GPU" = true ]; then
    echo "Skipping GPU (cupy) setup (--no-gpu given)."
else
    CUDA_MAJOR=""
    if command -v nvcc &> /dev/null; then
        CUDA_MAJOR="$(nvcc --version 2>/dev/null | grep -oP 'release \K[0-9]+' | head -1)"
    elif command -v nvidia-smi &> /dev/null; then
        # nvidia-smi reports the driver's supported CUDA version, not necessarily an installed
        # toolkit -- good enough to pick a cupy wheel, which bundles its own CUDA runtime libs.
        CUDA_MAJOR="$(nvidia-smi 2>/dev/null | grep -oP 'CUDA Version: \K[0-9]+' | head -1)"
    fi

    if [ -z "$CUDA_MAJOR" ]; then
        echo "No CUDA toolkit/driver detected (checked 'nvcc' and 'nvidia-smi') -- skipping cupy."
        echo "GPU acceleration is optional; both layers run fine on CPU. If you later add a GPU,"
        echo "re-run this script, or install cupy yourself: https://docs.cupy.dev/en/stable/install.html"
    else
        CUPY_PKG="cupy-cuda${CUDA_MAJOR}x"
        echo "Detected CUDA major version ${CUDA_MAJOR} -- installing ${CUPY_PKG} into the venv ..."
        "$VENV_PIP" install "$CUPY_PKG" || {
            echo "Warning: failed to install ${CUPY_PKG}. Check your CUDA version and install the"
            echo "matching build manually: https://docs.cupy.dev/en/stable/install.html"
        }
    fi
fi

echo
echo "Done. Activate the venv with:"
echo "  source $VENV_PATH/bin/activate"
echo "Then set USE_GPU=1 in your environment (or pass --gpu to reorganise_dir.py) to enable GPU acceleration."
