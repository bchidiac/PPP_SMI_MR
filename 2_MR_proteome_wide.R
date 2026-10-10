# ---- Step 2 | Proteome-wide cis-MR ------------------------------------------
# Two-sample cis-MR of UKB-PPP plasma pQTLs on 10 psychiatric outcomes
# (Methods; Supplementary Methods SM3, SM4):
#   - one instrument per conditionally independent cis signal, matched on the
#     allele pair UKB-PPP reported; signals shared by several proteins dropped
#   - LD proxies (r2 >= 0.8) from the same cis window for instruments absent
#     from an outcome
#   - Steiger filtering on the observed scale
#   - Wald ratio for one instrument, correlated IVW for two or more, with
#     Cochran's Q from the same model
#   - Bonferroni within outcome, one test per protein-outcome pair
#   - correlated MR-Egger, leave-one-out, pruned weighted median and
#     per-instrument Wald ratios for significant multi-instrument pairs
#
# Reads   1a exposure extract, 1b LD panel, 1c outcome GWAS, UKB-PPP cis signals
# Writes  MR tests, heterogeneity and sensitivity tables, harmonised
#         instruments, Manhattan and scatter plots, cross-outcome summary
# Feeds   3a, 3b, 4, 5, 7      Option  MR_F_MIN=10 for the F >= 10 analysis

print(version$version.string)

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(MendelianRandomization)
library(tidyverse)
library(ggrepel)
library(readxl)

# ---- Settings ----
env <- function(name, default) { v <- Sys.getenv(name); if (nzchar(v)) v else default }

data_dir      <- env("MR_DATA_DIR", "../MR_data")
gwas_dir      <- env("MR_GWAS_DIR", file.path(data_dir, "target_GWASes_updated"))
exposure_file <- env("MR_EXPOSURE_FILE", file.path(
  data_dir, "UKB_PPP_GWASes_European_Cis",
  "all_chromosomes_GRCh37_GRCh38_european_discovery.tsv.gz"))
signals_file  <- env("MR_SIGNALS_FILE", file.path(data_dir, "conditionally_independent_snps.xlsx"))
ld_panel      <- env("MR_LD_PANEL", file.path(data_dir, "LD_reference_panel_1000G_EUR_phase3/EUR"))
plink_bin     <- env("MR_PLINK", genetics.binaRies::get_plink_binary())   # PLINK 1.9
out_dir       <- env("MR_OUT_DIR", "results")
f_min         <- as.numeric(env("MR_F_MIN", "0"))
run_outcomes  <- env("MR_OUTCOMES", "")            # comma-separated display names; blank = all
run_date      <- env("MR_RUN_DATE", as.character(Sys.Date()))

alpha         <- 0.05
collinear_r2  <- 0.95  # instruments this correlated in the panel are treated as one signal
prune_r2      <- 0.1   # LD pruning threshold for the weighted median
set.seed(20260921)     # weighted-median bootstrap

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
scatter_root <- file.path(out_dir, "scatter_plots")
dir.create(scatter_root, showWarnings = FALSE)
log_file <- file.path(out_dir, paste0("run_log_", run_date, ".txt"))
logf <- function(...) {
  line <- paste0(...)
  message(line)
  cat(line, "\n", file = log_file, append = TRUE, sep = "")
}
logf("Run ", run_date, "; exposure ", basename(exposure_file), "; F >= ", f_min,
     "; LD panel ", ld_panel)

# ---- Helpers ----
# observed scale: r2 = z^2 / (z^2 + N)
steiger_keep <- function(dat) {
  r2 <- function(b, se, n) { z <- b / se; z^2 / (z^2 + n) }
  r2_exp <- r2(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)
  r2_out <- r2(dat$beta.outcome,  dat$se.outcome,  dat$samplesize.outcome)
  keep   <- is.na(r2_exp) | is.na(r2_out) | r2_exp > r2_out
  logf("  Steiger filtering removed ", sum(!keep), " of ", length(keep), " instruments")
  dat[keep, ]
}

# LD between instruments: panel variants matched by rsID or position and alleles
panel_bim <- fread(paste0(ld_panel, ".bim"), header = FALSE,
                   col.names = c("chr", "SNP", "cm", "pos", "A1", "A2"))
panel_bim[, `:=`(idx = .I, chr = sub("^chr", "", as.character(chr)))]
panel_n <- nrow(fread(paste0(ld_panel, ".fam"), header = FALSE))

match_panel <- function(snp, chr, pos, ea, oa) {
  q <- data.table(i = seq_along(snp), SNP = snp, chr = sub("^chr", "", as.character(chr)),
                  pos = as.integer(pos), ea = toupper(ea), oa = toupper(oa))
  same <- function(x) x[(ea == A1 & oa == A2) | (ea == A2 & oa == A1)]
  by_id  <- same(merge(q, panel_bim[, .(SNP, idx, A1, A2)], by = "SNP"))
  rest   <- q[!i %in% by_id$i]
  by_pos <- same(merge(rest, panel_bim[, .(chr, pos, idx, A1, A2)], by = c("chr", "pos")))
  hit <- rbind(by_id[, .(i, idx, A1, A2)], by_pos[, .(i, idx, A1, A2)])[!duplicated(i)]
  out <- q[, .(i, ea)][hit, on = "i"]
  res <- data.table(i = seq_along(snp), idx = NA_integer_, sign = NA_real_)
  res[out$i, `:=`(idx = out$idx, sign = ifelse(out$ea == out$A1, 1, -1))]
  res
}

read_genotypes <- function(idx) {
  nbytes <- ceiling(panel_n / 4)
  con <- file(paste0(ld_panel, ".bed"), "rb"); on.exit(close(con))
  stopifnot(identical(readBin(con, "raw", 3), as.raw(c(0x6c, 0x1b, 0x01))))
  decode <- c(2, NA, 1, 0)   # 00 hom A1, 01 missing, 10 het, 11 hom A2 -> A1 dosage
  G <- matrix(NA_real_, length(idx), panel_n)
  for (k in order(idx)) {
    seek(con, 3 + (idx[k] - 1) * nbytes)
    b <- as.integer(readBin(con, "raw", nbytes))
    codes <- as.vector(rbind(b %% 4, (b %/% 4) %% 4, (b %/% 16) %% 4, (b %/% 64) %% 4))[1:panel_n]
    G[k, ] <- decode[codes + 1]
  }
  G
}

wald_ratio <- function(d) {
  b <- d$beta.outcome / d$beta.exposure
  se <- abs(d$se.outcome / d$beta.exposure)          # first order, as TwoSampleMR
  list(b = b, se = se, pval = 2 * pnorm(-abs(b / se)))
}

mr_in <- function(d, R) mr_input(bx = d$beta.exposure, bxse = d$se.exposure,
                                 by = d$beta.outcome, byse = d$se.outcome,
                                 correlation = R, snps = d$SNP)

ivw_correlated <- function(d, R) {
  f <- mr_ivw(mr_in(d, R), model = "random", correl = TRUE)
  list(b = f@Estimate, se = f@StdError, pval = f@Pvalue,
       Q = f@Heter.Stat[1], Q_df = nrow(d) - 1, Q_pval = f@Heter.Stat[2])
}

ivw_uncorrelated <- function(d) {                   # standard IVW, for comparison
  w <- 1 / d$se.outcome^2
  b <- sum(w * d$beta.exposure * d$beta.outcome) / sum(w * d$beta.exposure^2)
  se_f <- sqrt(1 / sum(w * d$beta.exposure^2))
  rse <- sqrt(sum(w * (d$beta.outcome - b * d$beta.exposure)^2) / (nrow(d) - 1))
  se <- se_f * max(1, rse)
  list(b = b, se = se, pval = 2 * pnorm(-abs(b / se)))
}

# keep the stronger instrument of any pair with r2 >= r2_max
prune_ld <- function(d, R, r2_max) {
  ord <- order(-abs(d$beta.exposure / d$se.exposure))
  keep <- integer(0)
  for (j in ord) if (!length(keep) || all(R[j, keep]^2 < r2_max)) keep <- c(keep, j)
  sort(keep)
}

# a proxy must also be present in the outcome, and enters with its own
# exposure and outcome statistics
find_proxies <- function(missing_df, outcome_snps, r2_min = 0.8) {
  if (!nrow(missing_df)) return(tibble(gene = character(), SNP_A = character(),
                                       SNP_B = character(), r2 = numeric(), pval.exposure = numeric(), dist = numeric()))
  inst <- as.data.table(missing_df)[, .(gene, SNP_A = SNP, chr.exposure, pos_A = pos.exposure,
                                        effect_allele.exposure, other_allele.exposure)]
  m <- match_panel(inst$SNP_A, inst$chr.exposure, inst$pos_A,
                   inst$effect_allele.exposure, inst$other_allele.exposure)
  inst[, panel_id := panel_bim$SNP[m$idx]]
  inst <- inst[!is.na(panel_id)]
  tmp <- tempfile("proxies_")
  writeLines(unique(inst$panel_id), paste0(tmp, ".snps"))
  system2(plink_bin, c("--bfile", shQuote(ld_panel), "--r2", "--ld-snp-list", shQuote(paste0(tmp, ".snps")),
                       "--ld-window-kb", "5000", "--ld-window", "5000", "--ld-window-r2", r2_min,
                       "--out", shQuote(tmp)), stdout = FALSE, stderr = FALSE)
  ld <- fread(paste0(tmp, ".ld"))[SNP_A != SNP_B & R2 >= r2_min, .(panel_id = SNP_A, SNP_B, R2)]
  ld <- merge(ld, panel_bim[, .(SNP_B = SNP, chr_B = chr, pos_B = pos, A1_B = A1, A2_B = A2)], by = "SNP_B")
  cand <- merge(inst[, .(gene, SNP_A, pos_A, panel_id)], ld, by = "panel_id", allow.cartesian = TRUE)
  cand[, dist := abs(pos_B - pos_A)]
  cis <- as.data.table(cis_SNPs)[gene %in% cand$gene,
           .(gene, SNP, chr = as.character(chr.exposure), pos = as.integer(pos.exposure),
             key = paste(pmin(effect_allele.exposure, other_allele.exposure),
                         pmax(effect_allele.exposure, other_allele.exposure)), pval.exposure)]
  cand[, key := paste(pmin(A1_B, A2_B), pmax(A1_B, A2_B))]
  by_id  <- merge(cand, cis[, .(gene, SNP_B = SNP, SNP, pval.exposure)], by = c("gene", "SNP_B"))
  rest   <- cand[!paste(gene, SNP_A, SNP_B) %in% by_id[, paste(gene, SNP_A, SNP_B)]]
  by_pos <- merge(rest, cis[, .(gene, chr_B = chr, pos_B = pos, key, SNP, pval.exposure)],
                  by = c("gene", "chr_B", "pos_B", "key"))
  out <- rbind(by_id[, .(gene, SNP_A, SNP_B = SNP, r2 = R2, pval.exposure, dist)],
               by_pos[, .(gene, SNP_A, SNP_B = SNP, r2 = R2, pval.exposure, dist)])
  as_tibble(unique(out[SNP_B %in% outcome_snps & SNP_B != SNP_A]))
}

make_manhattan <- function(df_method,
                           top_genes,
                           threshold,
                           title_main = "",
                           title_sub  = "") {
  don <- df_method %>%
    group_by(chr.outcome) %>%
    summarise(chr_len = max(pos.outcome), .groups = "drop") %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%   # genome length overflows integers
    select(-chr_len) %>%
    left_join(df_method, by = "chr.outcome") %>%
    arrange(chr.outcome, pos.outcome) %>%
    mutate(
      BPcum        = pos.outcome + tot,
      is_highlight = significant,
      is_annotate  = id.exposure %in% top_genes
    )

  axisdf <- don %>%
    group_by(chr.outcome) %>%
    summarise(center = mean(range(BPcum)), .groups = "drop")

  ymax <- max(-log10(c(don$pval, threshold)), na.rm = TRUE)
  pad  <- ymax * 0.05

  ggplot(don, aes(x = BPcum, y = -log10(pval))) +
    geom_point(aes(colour = as.factor(chr.outcome)), size = 1.3, alpha = 0.8) +
    scale_colour_manual(values = rep(c("grey", "skyblue"), 22)) +
    geom_point(data = subset(don, is_highlight), colour = "orange", size = 1.3) +
    ggrepel::geom_label_repel(
      data = subset(don, is_annotate) %>% distinct(id.exposure, .keep_all = TRUE),
      aes(label = id.exposure),
      size = 6, max.overlaps = 20,
      box.padding = .35, point.padding = .5, segment.color = "grey50") +
    scale_x_continuous(breaks = axisdf$center, labels = axisdf$chr.outcome) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax + pad)) +
    geom_hline(yintercept = -log10(threshold), linetype = "dashed") +
    labs(x = "Chromosome",
         y = expression(-log[10](italic(p))),
         title = title_main, subtitle = title_sub) +
    theme_bw() +
    theme(
      legend.position    = "none",
      panel.border       = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x       = element_text(size = 14),
      axis.title.y       = element_text(size = 14),
      plot.title         = element_text(size = 20, hjust = .5),
      plot.subtitle      = element_text(size = 12, hjust = .5)
    )
}

# ---- 1. Instruments ----
# autosomes only: chrX genotype coding differs between GWAS
cis_SNPs <- data.table::fread(exposure_file) %>%
  filter(as.character(chr.exposure) %in% as.character(1:22)) %>%
  mutate(exposure = gene, SNP = id.exposure) %>%
  mutate(id.exposure = gene) %>%
  select(-pos.exposure.38) %>%
  mutate(pos.exposure = pos.exposure.19) %>%
  select(-pos.exposure.19)

dt <- read_excel(signals_file, skip = 4)
names(dt)[1] <- "Protein ID"
cis_signals <- dt %>%
  filter(`Cis/trans` == "cis") %>%
  mutate(gene = sub(":.*", "", `Protein ID`),
         allele_key = {
           parts <- strsplit(sub("^[^:]+:[^:]+:", "", `Variant ID`), ":")
           vapply(parts, function(a) paste(sort(a[1:2]), collapse = "/"), character(1))
         })

# at multi-allelic sites keep the allele pair UKB-PPP reported, or the rsID
# match where none does
filtered_cis_SNPs <- cis_SNPs %>%
  filter(paste0(gene, SNP) %in% paste0(cis_signals$gene, cis_signals$rsID)) %>%
  mutate(allele_key = paste(pmin(effect_allele.exposure, other_allele.exposure),
                            pmax(effect_allele.exposure, other_allele.exposure), sep = "/"),
         st16_match = paste0(gene, SNP, allele_key) %in%
                      paste0(cis_signals$gene, cis_signals$rsID, cis_signals$allele_key)) %>%
  group_by(gene, SNP) %>%
  filter(!any(st16_match) | st16_match) %>%
  ungroup() %>%
  select(-allele_key, -st16_match)

# proteins on several Olink panels: keep the largest N, then the smallest p
dedup_panels <- function(x) {
  x <- as.data.table(x)
  setorder(x, id.exposure, SNP, -samplesize.exposure, pval.exposure)
  as_tibble(unique(x, by = c("id.exposure", "SNP")))
}
filtered_cis_SNPs <- dedup_panels(filtered_cis_SNPs)

cis_SNPs <- dedup_panels(cis_SNPs)

# the second filter bars them from the proxy pool too
shared_cis_rsids <- cis_signals %>% distinct(gene, rsID) %>% count(rsID) %>%
  filter(n > 1) %>% pull(rsID)
n_before <- nrow(filtered_cis_SNPs)
filtered_cis_SNPs <- filtered_cis_SNPs %>% filter(!SNP %in% shared_cis_rsids)
cis_SNPs <- cis_SNPs %>% filter(!SNP %in% shared_cis_rsids)
logf("Instruments: ", n_before, " matched to cis signals; ", n_before - nrow(filtered_cis_SNPs),
     " removed as cis signals for more than one protein (", length(shared_cis_rsids),
     " shared variants); ", nrow(filtered_cis_SNPs), " kept for ",
     n_distinct(filtered_cis_SNPs$id.exposure), " proteins")

# ---- 2. MR per outcome ----
outcome_files <- c(
  "202205_updated_MVP_MGB_VU_normmgh.meta.0.8n.sumstat.tsv.gz",
  "OUDMTAG_Sumstats_Deak2022_MTAGcolumns_EURfreq_clean.tsv.gz",
  "AntiDepPercImprovEUR2021_phenotypeN.tsv.gz",
  "summary_stats_release_finngen_R12_F5_PANIC_GRCh37_rsID_split.tsv.gz",
  "bip2024_eur_noUKB_no23andMe.tsv.gz",
  "PGC3_SCZ_wave3_european_autosome_clean_MAF_GRCh37.tsv.gz",
  "pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv.gz",
  "daner_OCDmeta_wo23andMe_LOOUKBB_080425.tsv.gz",
  "eur_ptsd_pcs_v4_aug3_2021.tsv.gz",
  "trs_interaction_effectiveN.tsv.gz"
)
outcome_names <- c(
  "Treatment resistance in depression",
  "Opioid use disorder",
  "Antidepressant response",
  "Panic disorder",
  "Bipolar disorder",
  "Schizophrenia",
  "Major depressive disorder",
  "Obsessive compulsive disorder",
  "Post-traumatic stress disorder",
  "Treatment resistance in schizophrenia"
)
if (nzchar(run_outcomes)) {
  sel <- outcome_names %in% trimws(strsplit(run_outcomes, ",")[[1]])
  outcome_files <- outcome_files[sel]; outcome_names <- outcome_names[sel]
}

for (i in seq_along(outcome_files)) {
  outcome_file  <- outcome_files[i]
  outcome_name  <- outcome_names[i]
  outcome_short <- gsub("[.].*", "", outcome_file)
  stamp <- function(prefix, ext = "csv")
    file.path(out_dir, paste0(prefix, "_", outcome_short, "_", run_date, ".", ext))

  logf("\nOutcome: ", outcome_name, " (", outcome_file, ")")

  scatter_dir <- file.path(scatter_root, outcome_short)
  dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)

  target_GWAS <- data.table::fread(file.path(gwas_dir, outcome_file)) %>%
    mutate(id.outcome = outcome_name, outcome = outcome_name)

  non_clump_dat <- harmonise_data(exposure_dat = filtered_cis_SNPs,
                                  outcome_dat  = target_GWAS, action = 2)

  missing_df     <- filtered_cis_SNPs %>% filter(!SNP %in% non_clump_dat$SNP)
  unique_missing <- unique(missing_df$SNP)

  proxies <- find_proxies(missing_df, target_GWAS$SNP)

  # best proxy per instrument: highest r2, then smallest exposure p, then closest
  proxies <- proxies %>%
    group_by(gene, SNP_A) %>%
    arrange(desc(r2), pval.exposure, dist, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()

  proxies_SNPs <- cis_SNPs %>%
    semi_join(proxies %>% transmute(gene, SNP = SNP_B), by = c("gene", "SNP"))

  filtered_cis_SNPs_add_proxies <- bind_rows(filtered_cis_SNPs, proxies_SNPs)
  non_clump_dat <- harmonise_data(exposure_dat = filtered_cis_SNPs_add_proxies,
                                  outcome_dat  = target_GWAS, action = 2)
  logf("  Instruments missing from outcome: ", length(unique_missing),
       "; proxies found: ", nrow(proxies_SNPs))

  # drop ambiguous palindromes and duplicate rows, then Steiger (and F) filters
  dat <- as.data.table(non_clump_dat)[mr_keep == TRUE]
  n_dup <- sum(duplicated(dat, by = c("id.exposure", "SNP")))
  dat <- unique(dat, by = c("id.exposure", "SNP"))
  logf("  Harmonised instruments: ", nrow(dat), " (", n_dup, " duplicate protein-SNP rows removed)")
  dat <- as.data.table(steiger_keep(as.data.frame(dat)))
  dat[, F_stat := (beta.exposure / se.exposure)^2]
  if (f_min > 0) {
    logf("  F >= ", f_min, " filter removed ", sum(dat$F_stat < f_min), " instruments")
    dat <- dat[F_stat >= f_min]
  }

  m <- match_panel(dat$SNP, dat$chr.exposure, dat$pos.exposure,
                   dat$effect_allele.exposure, dat$other_allele.exposure)
  dat[, `:=`(panel_idx = m$idx, panel_sign = m$sign)]
  multi <- dat[, .N, by = id.exposure][N >= 2, id.exposure]
  need <- unique(dat[id.exposure %in% multi & !is.na(panel_idx), panel_idx])
  G <- read_genotypes(need); rownames(G) <- need
  # a variant with no variation in the panel cannot contribute a correlation, so it is treated as absent
  mono <- rownames(G)[apply(G, 1, function(g) sd(g, na.rm = TRUE)) %in% c(0, NA)]
  dat[as.character(panel_idx) %in% mono, panel_idx := NA]
  logf("  Instruments of multi-instrument proteins not in the LD panel: ",
       sum(dat$id.exposure %in% multi & is.na(dat$panel_idx)), " of ",
       sum(dat$id.exposure %in% multi))

  # instruments missing from the LD panel are left out, and a protein with none
  # of them in the panel falls back to its strongest; pairs with
  # r2 >= collinear_r2 count as one signal
  tests <- list(); kept_dat <- list()
  for (g in unique(dat$id.exposure)) {
    d <- dat[id.exposure == g]
    n_input <- nrow(d); n_absent <- 0L; n_collinear <- 0L; max_r2 <- NA_real_
    R <- NULL
    if (nrow(d) >= 2) {
      in_panel <- !is.na(d$panel_idx)
      n_absent <- sum(!in_panel)
      d <- if (any(in_panel)) d[in_panel] else d[which.max(abs(beta.exposure / se.exposure))]
      if (nrow(d) >= 2) {
        Gd <- G[as.character(d$panel_idx), , drop = FALSE]
        R <- suppressWarnings(cor(t(Gd), use = "pairwise.complete.obs")) * outer(d$panel_sign, d$panel_sign)
        R[is.na(R)] <- 0; diag(R) <- 1
        keep <- prune_ld(d, R, collinear_r2)
        n_collinear <- nrow(d) - length(keep)
        d <- d[keep]; R <- R[keep, keep, drop = FALSE]
        if (nrow(d) >= 2) max_r2 <- max(R[upper.tri(R)]^2)
      }
    }
    row <- list(id.exposure = g, id.outcome = outcome_name, outcome = outcome_name,
                nsnp_input = n_input, n_absent_panel = n_absent, n_collinear = n_collinear,
                nsnp = nrow(d), max_r2 = max_r2, min_F = min(d$F_stat))
    if (nrow(d) == 1) {
      e <- wald_ratio(d)
      row <- c(row, list(method = "Wald ratio", b = e$b, se = e$se, pval = e$pval,
                         Q = NA, Q_df = NA, Q_pval = NA,
                         b_ivw_uncorrelated = NA, se_ivw_uncorrelated = NA, pval_ivw_uncorrelated = NA))
    } else {
      e <- ivw_correlated(d, R); u <- ivw_uncorrelated(d)
      row <- c(row, list(method = "Inverse variance weighted", b = e$b, se = e$se, pval = e$pval,
                         Q = e$Q, Q_df = e$Q_df, Q_pval = e$Q_pval,
                         b_ivw_uncorrelated = u$b, se_ivw_uncorrelated = u$se,
                         pval_ivw_uncorrelated = u$pval))
    }
    tests[[g]] <- row
    kept_dat[[g]] <- list(d = d, R = R)
  }
  res <- rbindlist(tests)
  res[, exposure := id.exposure]   # column name read by the downstream scripts

  res[, n_tests := .N]
  res[, p_bonferroni := pmin(1, pval * n_tests)]
  res[, significant := p_bonferroni < alpha]
  logf("  Tests: ", nrow(res), " (", sum(res$method == "Wald ratio"), " Wald ratio, ",
       sum(res$method != "Wald ratio"), " IVW); instruments dropped as absent from panel: ",
       sum(res$n_absent_panel), "; as collinear: ", sum(res$n_collinear),
       "; Bonferroni threshold p < ", signif(alpha / nrow(res), 3), "; significant: ", sum(res$significant))

  # sensitivity estimators for significant multi-instrument associations
  sens <- list(); loo <- list(); two_snp <- list()
  for (g in res[significant == TRUE & nsnp >= 2, id.exposure]) {
    d <- kept_dat[[g]]$d; R <- kept_dat[[g]]$R
    if (nrow(d) == 2) {
      w <- wald_ratio(d)
      two_snp[[g]] <- data.table(id.exposure = g, id.outcome = outcome_name, SNP = d$SNP,
                                 b = w$b, se = w$se, pval = w$pval)
      next
    }
    eg <- mr_egger(mr_in(d, R), correl = TRUE)
    kp <- prune_ld(d, R, prune_r2)
    wm <- if (length(kp) >= 3)
      mr_median(mr_in(d[kp], R[kp, kp]), weighting = "weighted", iterations = 10000) else NULL
    lo <- rbindlist(lapply(seq_len(nrow(d)), function(j) {
      e <- ivw_correlated(d[-j], R[-j, -j, drop = FALSE])
      data.table(id.exposure = g, id.outcome = outcome_name, SNP_left_out = d$SNP[j],
                 b = e$b, se = e$se, pval = e$pval)
    }))
    loo[[g]] <- lo
    sens[[g]] <- data.table(
      id.exposure = g, id.outcome = outcome_name, nsnp = nrow(d),
      egger_b = eg@Estimate, egger_se = eg@StdError.Est, egger_pval = eg@Pvalue.Est,
      egger_intercept = eg@Intercept, egger_intercept_se = eg@StdError.Int,
      egger_intercept_pval = eg@Pvalue.Int,
      wm_nsnp = length(kp),
      wm_b = if (is.null(wm)) NA_real_ else wm@Estimate,
      wm_se = if (is.null(wm)) NA_real_ else wm@StdError,
      wm_pval = if (is.null(wm)) NA_real_ else wm@Pvalue,
      loo_max_pval = max(lo$pval), loo_same_sign = all(sign(lo$b) == sign(res[id.exposure == g, b])))
  }
  fwrite(rbindlist(sens), stamp("sensitivity"))
  fwrite(rbindlist(loo), stamp("leave_one_out"))
  fwrite(rbindlist(two_snp), stamp("two_instrument_wald_ratios"))

  fwrite(res, stamp("mr_results_all"))
  saveRDS(list(dat = dat, kept = kept_dat), stamp("harmonised", "rds"))
  fwrite(res[method != "Wald ratio", .(id.exposure, id.outcome, outcome, method, nsnp, Q, Q_df, Q_pval)],
         stamp("heterogeneity_stats"))

  res2 <- merge(res, dat[, .(gene, protein, SNP, beta.exposure, pval.exposure, beta.outcome,
                              pval.outcome, chr.outcome, pos.outcome)],
                by.x = "id.exposure", by.y = "gene", allow.cartesian = TRUE)
  fwrite(res2, stamp("plot_raw"))

  for (mth in c("Wald ratio", "Inverse variance weighted"))
    fwrite(res2[method == mth & significant == TRUE],
           stamp(paste0("text_MR_top_Gene_", gsub(" ", "_", mth))))

  top_genes <- res[significant == TRUE][order(pval)][1:min(.N, 20), id.exposure]
  plot <- make_manhattan(df_method = res2, top_genes = top_genes,
                         threshold = alpha / nrow(res),
                         title_main = sprintf("MR of UKB-PPP cis-pQTLs on\n%s", outcome_name),
                         title_sub = "Significant proteins (Bonferroni) in orange; top 20 labelled")
  ggsave(filename = stamp("Manhattan", "pdf"), plot = plot, device = "pdf", width = 10, height = 5)

  for (g in res[significant == TRUE & method != "Wald ratio", id.exposure]) {
    d <- as.data.frame(kept_dat[[g]]$d)
    p <- mr_scatter_plot(as.data.frame(res[id.exposure == g,
                           .(id.exposure, id.outcome, exposure = id.exposure, outcome, method, nsnp, b, se, pval)]),
                         d)[[1]]
    ggsave(filename = file.path(scatter_dir, sprintf("scatter_IVW_%s_%s_%s.pdf", g, outcome_short, run_date)),
           plot = p, width = 7, height = 5, device = "pdf")
  }
}

# ---- 3. Cross-outcome summary ----
all_res <- rbindlist(lapply(list.files(out_dir, paste0("^mr_results_all_.*_", run_date, "\\.csv$"),
                                       full.names = TRUE), fread))
summary_df <- all_res[significant == TRUE, .(
  n_disorders = uniqueN(id.outcome),
  disorders   = paste(sort(unique(id.outcome)), collapse = "; "),
  methods     = paste(unique(method), collapse = "; ")), by = .(protein = id.exposure)][order(-n_disorders)]
fwrite(summary_df, file.path(out_dir, paste0("proteins_significant_across_disorders_", run_date, ".csv")))
logf("\nSignificant protein-outcome pairs: ", sum(all_res$significant), "; proteins significant in >= 2 outcomes: ",
     sum(summary_df$n_disorders >= 2), "; in >= 3: ", sum(summary_df$n_disorders >= 3))
