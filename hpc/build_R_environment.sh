#!/bin/bash
#SBATCH -J build_coloc_env
#SBATCH -c 8
#SBATCH -t 03:00:00
# Conda environment used for the colocalisation runs: R 4.5.1, coloc 5.2.3, susieR 0.14.2,
# mixsqp 0.3.54, Rfast 2.1.5.2 and PLINK 2. Set ENV to the install location.
# Assumes mambaforge in the home directory; change the conda.sh path if yours differs.
set -euo pipefail
source ~/mambaforge/etc/profile.d/conda.sh
ENV=${ENV:-$PWD/envs/coloc_r451}
mamba create -y -p $ENV -c conda-forge -c bioconda r-base=4.5.1 plink2 \
  r-data.table r-dplyr r-tidyr r-ggplot2 r-ggrepel r-patchwork r-viridis r-readxl r-remotes r-r.utils \
  r-rcpp r-rcpparmadillo r-matrix r-irlba r-reshape r-crayon r-matrixstats r-rcppgsl gsl compilers make
conda activate $ENV
Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
remotes::install_version("mixsqp", "0.3.54", upgrade = "never")
remotes::install_version("susieR", "0.14.2", upgrade = "never")
remotes::install_version("Rfast", "2.1.5.2", upgrade = "never")
remotes::install_version("coloc", "5.2.3", upgrade = "never")
for (p in c("coloc", "susieR", "mixsqp", "Rfast", "data.table", "dplyr")) cat(p, as.character(packageVersion(p)), "\n")
cat(R.version.string, "\n")'
plink2 --version
conda env export -p $ENV > coloc_r451_environment.yml
Rscript -e 'writeLines(capture.output(sessionInfo()), "coloc_r451_R_session.txt")'
echo ENV_BUILD_DONE
