# ---- Step 0 | R packages ----------------------------------------------------
# Installs the R packages the pipeline uses. Run once.
#
# Also needs  PLINK 2 and bcftools on the PATH (steps 1b, 1c and 3a); Python 3
#             with pandas and synapseclient for step 1a. PLINK 1.9 (steps 1c
#             and 2) comes with genetics.binaRies.

cran <- c("data.table", "tidyverse", "dplyr", "tidyr", "readxl", "knitr", "rmarkdown",
          "ggplot2", "ggrepel", "ggtext", "viridis", "scales", "patchwork",
          "coloc", "susieR", "MendelianRandomization", "arrow", "openxlsx", "zip",
          "remotes", "BiocManager")
new <- setdiff(cran, rownames(installed.packages()))
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")

bioc <- c("rtracklayer", "GenomicRanges")
new <- setdiff(bioc, rownames(installed.packages()))
if (length(new)) BiocManager::install(new, update = FALSE, ask = FALSE)

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) remotes::install_github("MRCIEU/TwoSampleMR")
if (!requireNamespace("genetics.binaRies", quietly = TRUE)) remotes::install_github("MRCIEU/genetics.binaRies")
