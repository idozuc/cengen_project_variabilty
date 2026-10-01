#!/bin/bash
#SBATCH --job-name=basics_hvg
#SBATCH --time=16:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Submit with: sbatch --array=1-N --export=ALL,BASICS_RUN=/path/to/run slurm/run_basics_array.sh
set -euo pipefail

if [[ -z "${BASICS_RUN:-}" || -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "BASICS_RUN and SLURM_ARRAY_TASK_ID are required" >&2
  exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
srun Rscript "${PROJECT_ROOT}/scripts/basics/02_fit_cell_type.R" \
  --run-dir "${BASICS_RUN}" --task-id "${SLURM_ARRAY_TASK_ID}"
