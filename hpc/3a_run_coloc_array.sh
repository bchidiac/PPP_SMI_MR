#!/bin/bash
#SBATCH -J coloc
#SBATCH -c 2
#SBATCH -t 12:00:00
# Step 3a on a SLURM cluster, splitting the MR-significant pairs across tasks:
#   sbatch --array=0-<NTASKS-1> -o logs/coloc_%a.log 3a_run_coloc_array.sh settings.env
# settings.env sets the variables read by 3a_colocalisation.R
# (MR_DATA_DIR, MR_OUT_DIR, MR_RUN_DATE, COLOC_DIR, LD_MAT_DIR, PAIRS_FILE, NTASKS, ...).
# Activate an R environment with the packages from build_R_environment.sh first.
set -euo pipefail
set -a; source "$1"; set +a
export TMPDIR="${TMPDIR:-/tmp}/coloc_task_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
Rscript "${SCRIPTS:-.}/3a_colocalisation.R"
