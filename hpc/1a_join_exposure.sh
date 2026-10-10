#!/bin/bash
#SBATCH -J ukbppp_join
#SBATCH -c 8
#SBATCH -t 03:00:00
#SBATCH -o logs/join_%j.log
# Step 1a, second part: run after all extraction tasks have finished.
set -euo pipefail
OUT=${OUT:-ukbppp_european_extract}
SCRIPTS=${SCRIPTS:-.}
RSID_MAPS=${RSID_MAPS:-../MR_data/ukbppp_rsid_maps}
python "$SCRIPTS/1a_join_exposure_rsid.py" --indir "$OUT/cis_clean" --rsid-maps "$RSID_MAPS" \
  --out "$OUT/all_chromosomes_GRCh37_GRCh38_european_discovery.tsv.gz"
