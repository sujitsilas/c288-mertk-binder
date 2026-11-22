#!/usr/bin/env bash
set -euo pipefail
project=${1:?remote project directory}
campaign=${2:?remote campaign directory}
module purge
export PYTHONNOUSERSITE=1
cd "$project"
exec "$project/envs/project/bin/python" scripts/hpc/prepare_pilot_assessment.py "$campaign"
