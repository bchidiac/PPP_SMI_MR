# Proteome-wide *cis*-MR of plasma proteins across psychiatric outcomes

Code for the manuscript:

> **Functional genomic prioritisation of 1,788 plasma proteins across multiple mental illnesses and treatment response phenotypes.**

A proteome-wide two-sample *cis*-Mendelian randomisation (MR) study of UK Biobank Pharma Proteomics Project (UKB-PPP) plasma protein QTLs against 10 psychiatric outcomes (seven disorders and three treatment response phenotypes), with colocalisation of every MR-significant association.

## Pipeline

| Step | Script | Purpose |
|---|---|---|
| 0 | `install_dependencies.R` | Installs the R packages. |
| 1a | `1a_extract_exposure_cis.py`, `1a_join_exposure_rsid.py` | Downloads the UKB-PPP European discovery protein GWAS from Synapse, keeps *cis* variants (±1 Mb, INFO > 0.8, MAF > 0.05) and adds rsIDs and GRCh37/GRCh38 positions from the UKB-PPP rsID map. |
| 1b | `1b_build_LD_reference_panel.sh` | Builds the LD reference panel used throughout: the 503 European individuals of 1000 Genomes phase 3, autosomes, biallelic variants with MAF > 0.005, one record per variant ID (10,868,521 variants). |
| 1c | `1c_clean_outcome_GWAS.Rmd` | Cleans the ten outcome GWAS as downloaded (column names, derived fields, QC, and a liftover to GRCh37 for panic disorder) and counts their genome-wide significant loci (Table 1). |
| 2 | `2_MR_proteome_wide.R` | Main analysis: instruments, LD proxies, Steiger filtering, Wald ratio/correlated IVW, Bonferroni correction within outcome, sensitivity estimators, plots and a cross-outcome summary. |
| 3a | `3a_colocalisation.R` | `coloc.susie` for each MR-significant pair (*p*12 = 5×10⁻⁶ and 1×10⁻⁵), with `coloc.abf` where SuSiE finds no credible set. |
| 3b | `3b_coloc_instrument_overlap.R` | Whether MR instruments lie in (or in LD with) each colocalised 95% credible set. |
| 4 | `4_MR_method_and_heterogeneity_summary.R` | Estimator usage and heterogeneity across outcomes. |
| 5 | `5_cross_disorder_overlap.R` | Permutation test of whether proteins reach significance for several outcomes more often than chance. |
| 6 | `6_druggability.R` | Open Targets tractability and drugs, Human Protein Atlas secretome location, and a DGIdb cross-check. |
| 7 | `7_tables_and_figures.R` | Builds the manuscript tables and figures and the supplementary tables from the results. |

`hpc/` holds SLURM wrappers for the steps that are easier on a cluster (1a, 3a) and the R environment used for colocalisation.

## Running the pipeline

Steps 1a to 1c prepare the inputs and are run once:

```sh
# 1a  UKB-PPP cis extract, then the rsID join
python 1a_extract_exposure_cis.py --release discovery --outdir ukbppp_european
python 1a_join_exposure_rsid.py --indir ukbppp_european/cis_clean \
  --out ../MR_data/UKB_PPP_GWASes_European_Cis/all_chromosomes_GRCh37_GRCh38_european_discovery.tsv.gz

# 1b  LD reference panel (the script header gives the SLURM array form)
for c in $(seq 1 22); do SLURM_ARRAY_TASK_ID=$c bash 1b_build_LD_reference_panel.sh chr; done
bash 1b_build_LD_reference_panel.sh merge

# 1c  outcome GWAS cleaning, and the locus counts behind Table 1
Rscript -e 'rmarkdown::render("1c_clean_outcome_GWAS.Rmd")'
```

Steps 2 to 7 then run in order:

```sh
export MR_RUN_DATE=$(date +%F)         # stamps outputs, and selects the run each step reads
Rscript 2_MR_proteome_wide.R
Rscript 3a_colocalisation.R
Rscript 3b_coloc_instrument_overlap.R
Rscript 4_MR_method_and_heterogeneity_summary.R
Rscript 5_cross_disorder_overlap.R
Rscript 6_druggability.R
Rscript 7_tables_and_figures.R
```

Steps 2 and 4 to 6 write to `results/`, steps 3a and 3b to `coloc_results/`, and step 7 to `tables_and_figures/`.

The *F* ≥ 10 sensitivity analysis is step 2 run again with `MR_F_MIN=10` and its own `MR_OUT_DIR`; step 7 reads it from `MR_F10_DIR`.

Each script takes its paths from the settings block at its top, or from the environment variables named there, and otherwise expects the layout below:

| Variable | Default | Steps |
|---|---|---|
| `MR_RUN_DATE` | the current day; step 3a requires it | 2, 3a, 3b, 4, 5, 6, 7 |
| `MR_DATA_DIR` | `../MR_data` | 2, 3a, 6 |
| `MR_OUT_DIR` | `results` | 2, 3a, 3b, 4, 5, 6, 7 |
| `COLOC_DIR` | `coloc_results` | 3a, 3b, 7 |
| `MR_GWAS_DIR` | `../MR_data/target_GWASes_updated` | 2, 3a, 7 |
| `MR_LD_PANEL` | `../MR_data/LD_reference_panel_1000G_EUR_phase3/EUR` | 2, 3a, 3b |
| `MR_F_MIN` / `MR_F10_DIR` | `0` / `results_F10` | 2 / 7 |

## Data layout

The scripts expect a data folder next to the repository:

```
../MR_data/
  UKB_PPP_GWASes_European_Cis/          exposure file from step 1a
  ukbppp_rsid_maps/                     rsids_all_chromosomes_uniq2.tsv.gz (UKB-PPP rsID map, one
                                        record per variant ID: ID, rsid, POS19, POS38)
  conditionally_independent_snps.xlsx   UKB-PPP conditionally independent signals
  olink_protein_map_3k_v1.tsv           Olink protein-to-gene map (UKB-PPP release)
  raw_outcome_GWAS/                     outcome GWAS as downloaded
  target_GWASes_updated/                cleaned outcome GWAS from step 1c
  LD_reference_panel_1000G_EUR_phase3/  EUR.bed/.bim/.fam from step 1b
  liftover_chain_file/                  hg38ToHg19.over.chain (UCSC)
  druggability/
    OpenTargets_26.06/                  Open Targets Platform release 26.06 (parquet)
    HPA/                                proteinatlas.tsv.zip (Human Protein Atlas, accessed 21 September 2026)
    DGIdb/                              interactions.tsv (DGIdb 5.0)
```

## Requirements

R ≥ 4.3 (analyses were run in R 4.5.1), PLINK 2 and bcftools, and for step 1a Python 3 with `pandas` and `synapseclient` (a Synapse account is needed for the UKB-PPP data). `install_dependencies.R` installs the R packages; the exact versions used for colocalisation are pinned in `hpc/build_R_environment.sh`.

## Data availability

Summary statistics are not redistributed here.

- Exposure: UKB-PPP European (discovery) summary statistics (Sun et al. 2023), Synapse folder syn51365303: https://www.synapse.org/Synapse:syn51365303
- Outcomes: the GWAS listed in Table 1 of the manuscript, from their original sources.
- LD reference: built by step 1b from the 1000 Genomes phase 3 v5a VCFs (release 20130502), which carry rsIDs: https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/

## Licence

MIT; see [LICENSE](LICENSE).
