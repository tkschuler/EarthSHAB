#!/usr/bin/env bash
#
# Start gfs_collector.py in a detached, named `screen` session.
#
# Mirrors the "Running the collector in the background with screen" guide
# (docs/source/examples/downloadGFS.rst): it launches the long-running scheduled
# downloader inside a detached screen session named "gfs_collector" so it keeps
# polling NOAA's AWS GFS bucket after you log out.
#
# Usage:
#   scripts/start_gfs_collector.sh            # start the collector (no-op if already running)
#   CONDA_ENV=myenv scripts/start_gfs_collector.sh
#
# Then:
#   screen -r gfs_collector   # reattach to watch live output  (Ctrl-a d to detach)
#   screen -ls                # list sessions
#
set -euo pipefail

SESSION="gfs_collector"
CONDA_ENV="${CONDA_ENV:-EarthSHAB}"

# Project root = parent of this script's directory, so the script works no
# matter where it's invoked from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if ! command -v screen >/dev/null 2>&1; then
    echo "ERROR: 'screen' is not installed (sudo apt install screen / brew install screen)." >&2
    exit 1
fi

# Don't stack a second collector on top of a running one.
if screen -ls 2>/dev/null | grep -q "[.]${SESSION}[[:space:]]"; then
    echo "Session '${SESSION}' is already running. Reattach with: screen -r ${SESSION}"
    exit 0
fi

# Detached named session that activates the env, cd's to the repo, and runs the
# collector module — the one-liner form from the docs' tip.
#
# `conda activate` needs conda's shell function, which a non-interactive shell
# doesn't load — so we source conda.sh first. Without this the activate fails
# ("Run 'conda init' before 'conda activate'"), the && chain short-circuits, and
# the session dies before Python ever starts (with no visible error out here).
screen -dmS "$SESSION" bash -lc \
    "source \"\$(conda info --base)/etc/profile.d/conda.sh\" && conda activate '${CONDA_ENV}' && cd '${REPO_ROOT}' && python -m EarthSHAB.gfs_collector"

echo "Started '${SESSION}' (conda env: ${CONDA_ENV}, cwd: ${REPO_ROOT})."
echo "Reattach to watch output:  screen -r ${SESSION}   (detach with Ctrl-a d)"
echo "Stop it:                   screen -r ${SESSION}, then Ctrl-c, then exit"
