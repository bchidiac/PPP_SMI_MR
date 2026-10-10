#!/bin/bash
#SBATCH -J ukbppp_cis
#SBATCH -c 2
#SBATCH -t 06:00:00
#SBATCH --array=0-19
#SBATCH -o logs/extract_%A_%a.log
# Step 1a on a SLURM cluster: 20 tasks share the ~2,900 UKB-PPP protein archives.
# Set the account/partition for your cluster, and OUT and SCRIPTS below.
# Create the logs/ directory before submitting: SLURM resolves -o before the job starts.
# Python needs pandas and synapseclient.
set -euo pipefail
OUT=${OUT:-ukbppp_european_extract}
SCRIPTS=${SCRIPTS:-.}
python "$SCRIPTS/1a_extract_exposure_cis.py" --release discovery --outdir "$OUT" \
  --task "$SLURM_ARRAY_TASK_ID" --ntasks 20
