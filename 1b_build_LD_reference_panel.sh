#!/bin/bash
#SBATCH -J 1000g_eur
#SBATCH -c 4
#SBATCH -t 06:00:00
# ---- Step 1b | LD reference panel -------------------------------------------
# 1000 Genomes phase 3 LD reference, GRCh37 (Supplementary Methods SM4):
#   - the 503 European (EUR) individuals, autosomes
#   - biallelic records (SNVs, indels, structural variants), MAF > 0.005
#   - unique IDs: variants without an rsID get a chr:pos:ref:alt ID, and the
#     split records of multi-allelic sites, which share an rsID, are dropped
#
# Used by  1c, 2, 3a and 3b
# Needs    bcftools and plink2 on the PATH. Create logs/ before submitting:
#          SLURM resolves -o before the job starts.
#
#   sbatch --array=1-22 -o logs/chr_%a.log 1b_build_LD_reference_panel.sh chr
#   sbatch --dependency=afterok:<jobid> -o logs/merge.log 1b_build_LD_reference_panel.sh merge
#
# Without SLURM: for c in $(seq 1 22); do SLURM_ARRAY_TASK_ID=$c bash 1b_build_LD_reference_panel.sh chr; done
set -euo pipefail
OUT=${OUT:-../MR_data/LD_reference_panel_1000G_EUR_phase3}
URL=https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502   # v5a: same genotypes as EBI's v5b files, but with rsIDs
PANEL=integrated_call_samples_v3.20130502.ALL.panel
mkdir -p "$OUT/per_chr"; cd "$OUT"

if [ "$1" = "chr" ]; then
  c=$SLURM_ARRAY_TASK_ID
  VCF=ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  [ -s "$PANEL" ] || wget -q "$URL/$PANEL" || true
  awk '$3 == "EUR" {print $1}' "$PANEL" > eur_samples.txt
  { printf "#FID\tIID\tSEX\n"; awk '$3 == "EUR" {print "0\t" $1 "\t" ($4 == "male" ? 1 : 2)}' "$PANEL"; } > eur_sex.txt
  cd per_chr
  wget -q -c --tries=30 --waitretry=30 --timeout=120 "$URL/$VCF"
  gzip -t "$VCF"
  md5sum "$VCF" > "chr${c}.vcf.md5"   # records which download the panel was built from
  plink2 --vcf "$VCF" --keep ../eur_samples.txt --const-fid 0 --update-sex ../eur_sex.txt \
    --max-alleles 2 --maf 0.005 --set-missing-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 \
    --rm-dup exclude-all --make-bed --out "EUR_chr${c}"
  rm -f "$VCF"
elif [ "$1" = "merge" ]; then
  cd per_chr
  for c in $(seq 1 22); do echo "EUR_chr${c}"; done > merge_list.txt
  plink2 --pmerge-list merge_list.txt bfile --make-bed --out ../EUR
  cd .. && md5sum EUR.bed EUR.bim EUR.fam > EUR.md5
fi
