# ---- Step 3a | Colocalisation ----------------------------------------------
# coloc.susie for each MR-significant protein-outcome pair (Supplementary
# Methods SM5):
#   - an LD matrix per pair built with PLINK 2, with exposure and outcome
#     harmonised to the reference alleles
#   - run at each p12 in P12_GRID (5e-6 main, 1e-5 sensitivity)
#   - coloc.abf, which assumes a single causal variant, where SuSiE returns no
#     credible set or does not converge
#   - SuSiE run with an iteration cap of 10,000
#
# Reads   text_MR_top_Gene_* from step 2, 1a exposure extract, 1b LD panel
#         and 1c outcome GWAS
# Writes  per-pair summaries, per-SNP posteriors (both methods), regional plots
#         and a status table
# Feeds   3b, 7      Cluster  hpc/3a_run_coloc_array.sh splits pairs

library(data.table); library(dplyr); library(tidyr); library(coloc)
library(ggplot2);    library(ggrepel); library(viridis)
library(patchwork)

safe_dir <- function(x) gsub("\\s+", "_", x)

# ---- Settings ----
env <- function(x, default = NULL) {
  v <- Sys.getenv(x)
  if (nzchar(v)) return(v)
  if (is.null(default)) stop("missing setting ", x)
  default
}
data_dir      <- env("MR_DATA_DIR", "../MR_data")
exposure_file <- env("MR_EXPOSURE_FILE", file.path(data_dir, "UKB_PPP_GWASes_European_Cis",
                     "all_chromosomes_GRCh37_GRCh38_european_discovery.tsv.gz"))
mr_dir     <- env("MR_OUT_DIR", "results")
gwas_dir   <- env("MR_GWAS_DIR", file.path(data_dir, "target_GWASes_updated"))
ld_prefix  <- env("MR_LD_PANEL", file.path(data_dir, "LD_reference_panel_1000G_EUR_phase3/EUR"))
ld_mat_dir <- env("LD_MAT_DIR", "LD_matrices")
out_dir    <- env("COLOC_DIR", "coloc_results")
run_date   <- env("MR_RUN_DATE")
p12_grid   <- as.numeric(strsplit(env("P12_GRID", "5e-6,1e-5"), ",")[[1]])
susie_maxit <- 10000
task       <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")); ntasks <- as.integer(env("NTASKS", "1"))
pairs_file <- env("PAIRS_FILE", "")                     # optional: tab-separated outcome, gene
task_keys  <- NULL
if (nzchar(pairs_file)) {
  pairs     <- data.table::fread(pairs_file, sep = "\t")
  pairs     <- pairs[(seq_len(nrow(pairs)) - 1) %% ntasks == task]
  task_keys <- paste(pairs$outcome, pairs$gene)
}
# outcomes analysed as quantitative traits (sdY = 1 where betas are already standardised)
quant_sdY  <- list("Opioid use disorder" = 1, "Post-traumatic stress disorder" = 1,
                   "Treatment resistance in depression" = NA, "Antidepressant response" = NA)
status_file <- file.path(out_dir, sprintf("pair_status_task%03d.csv", task))
record <- function(outcome, gene, p12, status, abf = NULL) {
  row <- data.table::data.table(outcome = outcome, gene = gene, p12 = p12, status = status,
    abf_nsnps = if (is.null(abf)) NA else abf[["nsnps"]], abf_h3 = if (is.null(abf)) NA else abf[["PP.H3.abf"]],
    abf_h4 = if (is.null(abf)) NA else abf[["PP.H4.abf"]])
  data.table::fwrite(row, status_file, append = file.exists(status_file))
}
cat(sprintf("task %d/%d: %s pairs; p12 %s; susie maxit %d\n", task, ntasks,
            if (is.null(task_keys)) "all" else length(task_keys),
            paste(p12_grid, collapse = ","), susie_maxit))
thr_keep   <- 0.8                       # PPH4 threshold for colocalisation
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cis_SNPs <- fread(exposure_file) %>%
  mutate(
    exposure     = gene,
    SNP          = id.exposure,      # adds the SNP column the loop expects
    id.exposure  = gene,
    pos.exposure = pos.exposure.38   # kept for ordering and the regional plots
  ) %>%
  select(-pos.exposure.38, -pos.exposure.19)

# ---- Inputs ----
outcome_map <- c(
  "202205_updated_MVP_MGB_VU_normmgh"                   = "Treatment resistance in depression",
  "OUDMTAG_Sumstats_Deak2022_MTAGcolumns_EURfreq_clean" = "Opioid use disorder",
  "AntiDepPercImprovEUR2021_phenotypeN"                      = "Antidepressant response",
  "summary_stats_release_finngen_R12_F5_PANIC_GRCh37_rsID_split" = "Panic disorder",
  "bip2024_eur_noUKB_no23andMe"                         = "Bipolar disorder",
  "PGC3_SCZ_wave3_european_autosome_clean_MAF_GRCh37"   = "Schizophrenia",
  "pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11"        = "Major depressive disorder",
  "daner_OCDmeta_wo23andMe_LOOUKBB_080425"              = "Obsessive compulsive disorder",
  "eur_ptsd_pcs_v4_aug3_2021"                           = "Post-traumatic stress disorder",
  "trs_interaction_effectiveN"                                     = "Treatment resistance in schizophrenia"
)
gwas_file <- function(stub) {
  f <- file.path(gwas_dir, paste0(stub, ".tsv.gz"))
  if (!file.exists(f)) f <- file.path(gwas_dir, paste0(stub, ".sumstats.tsv.gz"))
  f
}

# ---- Helpers ----
read_mr <- function(path){
  if (!file.exists(path)) { warning("file not found: ", path); return(data.frame(exposure = character())) }
  dat <- fread(path)
  if (!"exposure" %in% names(dat)) { warning("no 'exposure' column in ", path); return(data.frame(exposure = character())) }
  dat$exposure <- as.character(dat$exposure)
  dat <- dat[!is.na(dat$exposure), ]
  if (nrow(dat) == 0) dat <- data.frame(exposure = character())
  dat
}

read_target_gwas <- function(path, outcome_label){
  if(!file.exists(path)) stop("GWAS file not found: ", path)
  gwas <- fread(path)
  rename_map <- c(chr.outcome = "chr", pos.outcome = "position", SNP = "rsid",
                  other_allele.outcome = "oa", effect_allele.outcome = "ea",
                  eaf.outcome = "eaf.outcome", samplesize.outcome = "N",
                  beta.outcome = "beta", se.outcome = "se", pval.outcome = "pval")
  keep <- intersect(names(rename_map), names(gwas))
  setnames(gwas, keep, rename_map[keep])
  gwas[, MAF := pmin(eaf.outcome, 1 - eaf.outcome)]
  mandatory <- c("chr","position","rsid","oa","ea","beta","se","pval","N","MAF")
  missing <- setdiff(mandatory, names(gwas))
  if(length(missing)) stop("GWAS file ", path," lacks: ", paste(missing, collapse=", "))
  gwas <- gwas[, ..mandatory]
  gwas[, outcome := outcome_label]
  gwas
}

make_assoc_panel <- function(dat, R, lead, chr, start, end,
                             p_col, pos_col, rsid_col = "rsid",
                             title = "",
                             beta_col = NULL, se_col = NULL,
                             y_buffer = 0.08, y_cap = NULL) {
  stopifnot(lead %in% rownames(R), lead %in% colnames(R))
  stopifnot(all(c(rsid_col, pos_col, p_col) %in% names(dat)))
  dat2 <- dat %>% dplyr::filter(.data[[pos_col]] >= start, .data[[pos_col]] <= end)
  if (!is.null(beta_col) && !is.null(se_col) && all(c(beta_col, se_col) %in% names(dat2))) {
    z <- abs(dat2[[beta_col]] / dat2[[se_col]])
    logp_two_sided <- pnorm(-z, log.p = TRUE) + log(2)   # finite for huge z
    dat2$logp <- -(logp_two_sided) / log(10)
  } else {
    p <- pmin(pmax(dat2[[p_col]], .Machine$double.xmin), 1)
    dat2$logp <- -log10(p)
  }
  r2 <- rep(NA_real_, nrow(dat2))
  inR <- dat2[[rsid_col]] %in% rownames(R)
  r2[inR] <- R[lead, dat2[[rsid_col]][inR]]^2
  df <- dat2 %>% dplyr::mutate(r2 = r2, lead = .data[[rsid_col]] == lead)
  ymax  <- if (is.null(y_cap)) max(df$logp, na.rm = TRUE) else y_cap
  y_top <- ymax * (1 + y_buffer)
  p <- ggplot(df, aes(x = .data[[pos_col]]/1e6, y = logp)) +
    geom_point(aes(fill = r2), shape = 21, size = 2, colour = "black", alpha = 0.85, na.rm = TRUE) +
    scale_fill_viridis_c(option = "plasma", limits = c(0,1), breaks = seq(0,1,0.2),
                         name = expression(r^2), na.value = "grey85",
                         guide = guide_colourbar(barheight = unit(3, "cm"))) +
    labs(x = paste0("Position on Chr", chr, " (Mb)"),
         y = expression(-log[10](p)), title = title) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18), plot.title = element_text(size = 18)) +
    coord_cartesian(ylim = c(0, y_top)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)))
  p <- p + ggrepel::geom_text_repel(
    data = subset(df, lead), aes(label = .data[[rsid_col]]),
    min.segment.length = 0, segment.color = "black", segment.size = 0.4,
    box.padding = 0.3, point.padding = 0.6, max.overlaps = Inf, size = 4)
  p
}

plot_region_triptych <- function(d_h, R, chr, lead, center_pos = NULL, flank = 250e3,
                                 rsid_col = "rsid", pos_col = "position",
                                 p_exp_col = "pval", p_out_col = "i.pval",
                                 exposure_label = "Exposure", outcome_label = "Outcome") {
  stopifnot(all(c(rsid_col, pos_col, p_exp_col, p_out_col) %in% names(d_h)))
  stopifnot(lead %in% d_h[[rsid_col]], lead %in% rownames(R))
  if (is.null(center_pos)) center_pos <- d_h[[pos_col]][match(lead, d_h[[rsid_col]])]
  start <- center_pos - flank; end <- center_pos + flank
  p1 <- make_assoc_panel(dat = d_h, R = R, lead = lead, chr = chr, start = start, end = end,
    p_col = p_exp_col, pos_col = pos_col, rsid_col = rsid_col, title = paste(exposure_label),
    beta_col = "beta", se_col = "se")
  p2 <- make_assoc_panel(dat = d_h, R = R, lead = lead, chr = chr, start = start, end = end,
    p_col = p_out_col, pos_col = pos_col, rsid_col = rsid_col, title = paste(outcome_label),
    beta_col = "i.beta", se_col = "i.se")
  (p1 / p2) + plot_layout(heights = c(3,3))
}

# coloc names this column SNP.PP.H4.row<n> when there are several signal pairs
# and SNP.PP.H4.abf when there is one
pick_pp_col <- function(res_dt, r) {
  cands <- c(paste0("SNP.PP.H4.row", r), paste0("SNP.PP.H4.row", r, ".abf"),
             "SNP.PP.H4.abf", "SNP.PP.H4")
  cands[match(cands, names(res_dt), nomatch = 0L) > 0][1]
}

as_coloc_summary_dt <- function(s) {
  if (is.null(s)) return(data.table())
  if (is.data.frame(s) || is.matrix(s)) return(as.data.table(s))
  if (is.atomic(s) && !is.null(names(s))) return(as.data.table(as.list(s)))
  if (is.list(s) && !is.null(names(s)))   return(as.data.table(s))
  data.table()
}

# ---- Colocalisation ----
pvar <- data.table::fread(paste0(ld_prefix, ".bim"),
  col.names = c("CHR","ID","CM","BP","A1","A2"))[, .(rsid = ID, REF = A1, ALT = A2)]

config <- tibble(
  stub    = names(outcome_map),
  outcome = unname(outcome_map),
  mr_ivw  = file.path(mr_dir, paste0("text_MR_top_Gene_Inverse_variance_weighted_", stub, "_", run_date, ".csv")),
  mr_wald = file.path(mr_dir, paste0("text_MR_top_Gene_Wald_ratio_", stub, "_", run_date, ".csv")),
  gwas    = vapply(stub, gwas_file, character(1))
)

for (k in seq_len(nrow(config))) {
  conf <- config[k, ]
  message("=== processing ", conf$outcome, " ===")

  mr_genes <- bind_rows(
    read_mr(conf$mr_ivw)  %>% select(exposure),
    read_mr(conf$mr_wald) %>% select(exposure)) %>% distinct()
  if (!is.null(task_keys)) mr_genes <- mr_genes %>% filter(paste(conf$outcome, exposure) %in% task_keys)
  if (nrow(mr_genes) == 0) next

  gwas <- read_target_gwas(conf$gwas, conf$outcome)

  dis_dir <- file.path(out_dir, safe_dir(conf$outcome))
  dir.create(dis_dir, showWarnings = FALSE, recursive = TRUE)

  for (j in seq_along(p12_grid)) {
    p12_val <- p12_grid[j]

    for (i in seq_len(nrow(mr_genes))) {
      g     <- mr_genes$exposure[i]
      expGW <- cis_SNPs %>% filter(gene == g) %>%
        mutate(MAF = pmin(eaf.exposure, 1 - eaf.exposure)) %>%
        distinct(SNP, .keep_all = TRUE)
      if (nrow(expGW) == 0) next

      names(expGW) <- c("chr","id.exposure","oa","ea","eaf.exposure","INFO","N",
                        "beta","se","pval","dataset","gene","protein","exposure",
                        "rsid","position","MAF")
      expGW <- expGW %>% select(chr, position, ea, oa, beta, se, pval, rsid, N, MAF)

      # exposure and outcome variants shared by rsID (at most ~10,000 around the gene)
      d <- expGW[gwas, on = "rsid", nomatch = NULL]
      d <- d[!is.na(MAF) & !is.na(i.MAF)]
      if (nrow(d) == 0) next
      if (nrow(d) > 9999) { o <- floor((nrow(d) - 9999)/2); d <- d[(o + 1):(o + 9999)] }
      d <- d %>% arrange(position)

      this_chr <- unique(expGW$chr)
      ld_dir    <- file.path(ld_mat_dir, safe_dir(conf$outcome))
      rsids_fn  <- file.path(ld_dir, paste0(g, "_chr", this_chr, ".unphased.vcor1.vars"))
      mat_fn    <- sub("\\.vars$", "", rsids_fn)

      if (!file.exists(rsids_fn) || !file.exists(mat_fn)) {
        if (!grepl("^([1-9]|1[0-9]|2[0-2])$", this_chr)) { record(conf$outcome, g, p12_val, "not_autosomal"); next }
        dir.create(ld_dir, showWarnings = FALSE, recursive = TRUE)
        d_list <- expGW[gwas, on = "rsid", nomatch = NULL][!is.na(MAF) & !is.na(i.MAF)]
        if (nrow(d_list) > 9999) d_list <- d_list[(floor((nrow(d_list) - 9999)/2) + 1):(floor((nrow(d_list) - 9999)/2) + 9999)]
        snp_file <- file.path(ld_dir, paste0("snps_", g, ".txt"))
        writeLines(d_list$rsid, snp_file)
        plink_log <- system2("plink2", c("--bfile", shQuote(ld_prefix), "--chr", this_chr, "--extract", shQuote(snp_file),
                            "--rm-dup", "force-first", "--r-unphased", "square", "ref-based",
                            "--out", shQuote(file.path(ld_dir, paste0(g, "_chr", this_chr)))),
                stdout = TRUE, stderr = TRUE)
        if (!file.exists(rsids_fn) || !file.exists(mat_fn)) message("plink2 failed for ", g, ": ", paste(tail(plink_log, 3), collapse = " | "))
        if (!file.exists(rsids_fn) || !file.exists(mat_fn)) { record(conf$outcome, g, p12_val, "ld_failed"); next }
      }

      d_h  <- merge(d, pvar, by = "rsid", all.x = TRUE)
      vars <- scan(rsids_fn, what = character(), quiet = TRUE)
      d_h  <- d_h %>%
        dplyr::filter(rsid %in% vars, !is.na(REF), !is.na(ALT), !is.na(ea), !is.na(i.ea))

      # keep variants whose allele pairs match the reference, then align to its first allele
      same_pair <- function(a, b) (a == d_h$REF & b == d_h$ALT) | (a == d_h$ALT & b == d_h$REF)
      n_before <- nrow(d_h)
      d_h <- d_h[same_pair(d_h$ea, d_h$oa) & same_pair(d_h$i.ea, d_h$i.oa), ]
      message(sprintf("  %s: %d variants matched to the reference; %d removed by the allele check", g, n_before, n_before - nrow(d_h)))
      flip_exp <- d_h$ea   != d_h$REF
      flip_out <- d_h$i.ea != d_h$REF
      d_h$beta[flip_exp]   <- -d_h$beta[flip_exp]
      d_h$i.beta[flip_out] <- -d_h$i.beta[flip_out]

      is_palin <- (d_h$REF %in% c("A","T")) & (d_h$ALT %in% c("A","T")) |
                  (d_h$REF %in% c("C","G")) & (d_h$ALT %in% c("C","G"))
      d_h <- d_h[!is_palin & d_h$MAF > 0.05, ]
      d_h <- d_h %>% arrange(position)

      ld_rsids  <- fread(rsids_fn, header = FALSE)
      ld_matrix <- fread(mat_fn)
      snps   <- ld_rsids[[1]]
      ld_mat <- as.matrix(ld_matrix)
      rownames(ld_mat) <- snps; colnames(ld_mat) <- snps

      keep_snps <- d_h$rsid
      ld_mat   <- ld_mat[keep_snps, keep_snps, drop = FALSE]
      ld_rsids <- ld_rsids[V1 %in% keep_snps]
      stopifnot(identical(rownames(ld_mat), keep_snps))

      rs_d  <- as.character(d_h$rsid)
      rs_ld <- as.character(ld_rsids$V1)
      if (length(rs_d) != length(rs_ld))
        stop(sprintf("Length mismatch: d has %d, ld has %d", length(rs_d), length(rs_ld)))

      # both traits take the outcome's GRCh37 positions; variants are aligned by rsID
      D1 <- list(type = "quant", beta = d_h$beta, varbeta = d_h$se^2, pvalues = d_h$pval,
                 snp = d_h$rsid, N = d_h$N[1], MAF = d_h$MAF, position = d_h$i.position, LD = ld_mat)
      D2 <- list(type = "cc", beta = d_h$i.beta, varbeta = d_h$i.se^2, pvalues = d_h$i.pval,
                 snp = d_h$rsid, N = d_h$i.N[1], MAF = d_h$i.MAF, position = d_h$i.position, LD = ld_mat)
      if (conf$outcome %in% names(quant_sdY)) {
        D2$type <- "quant"
        if (!is.na(quant_sdY[[conf$outcome]])) D2$sdY <- quant_sdY[[conf$outcome]]
      }
      run_abf <- function() {
        a <- coloc.abf(D1[setdiff(names(D1), "LD")], D2[setdiff(names(D2), "LD")], p12 = p12_val)
        data.table::fwrite(as.data.table(a$results)[, .(snp, pp = SNP.PP.H4)], file.path(dis_dir,
          sprintf("abf_coloc_snps_%s_%s_p12_%s_%s.csv", g, gsub("\\s+", "_", conf$outcome), format(p12_val), run_date)))
        a$summary
      }

      coloc_res <- tryCatch(
        coloc.susie(D1, D2, p12 = p12_val,
                    susie.args = list(maxit = susie_maxit, repeat_until_convergence = FALSE)),
        error = function(e) { message("coloc.susie failed for ", g, ": ", e$message); NULL })
      if (is.null(coloc_res)) { record(conf$outcome, g, p12_val, "susie_not_converged", run_abf()); next }

      summ <- as_coloc_summary_dt(coloc_res$summary)
      if (!nrow(summ) || !("PP.H4.abf" %in% names(summ))) {
        message("  -> no PP.H4.abf returned for ", g, " (skipping)")
        record(conf$outcome, g, p12_val, "no_credible_set", run_abf()); next
      }
      summ[, row := .I]
      record(conf$outcome, g, p12_val, "susie_result")

      summ_out <- copy(summ)[, `:=`(gene = g, outcome = conf$outcome,
                                    chr = this_chr, p12 = p12_val)]
      data.table::fwrite(summ_out, file = file.path(dis_dir,
        sprintf("susie_coloc_summary_%s_%s_p12_%s_%s.csv", g, gsub("\\s+","_",conf$outcome), format(p12_val), run_date)))

      res_dt <- as.data.table(coloc_res$results)
      for (r in seq_len(nrow(summ))) {
        ppcol <- pick_pp_col(res_dt, r)
        if (is.na(ppcol) || !nzchar(ppcol)) {
          message("  -> no per-SNP PP4 column for row ", r); next
        }
        rs <- res_dt[, .(snp, pp = get(ppcol))][is.finite(pp) & !is.na(pp)]
        pos_map <- unique(d_h[, .(snp = rsid, position = position)])
        rs <- left_join(rs, pos_map)
        rs_out <- rs[order(-pp)]
        data.table::fwrite(rs_out, file = file.path(dis_dir,
          sprintf("susie_coloc_snps_%s_%s_p12_%s_row%i_%s.csv",
                  g, gsub("\\s+","_",conf$outcome), format(p12_val), r, run_date)))

        if (isTRUE(summ$PP.H4.abf[r] >= thr_keep)) {
          lead_shared <- rs_out$snp[1]
          center_pos <- rs_out$position[1]
          if (is.na(center_pos)) center_pos <- NULL
          plt <- plot_region_triptych(
            d_h, ld_mat, chr = this_chr, lead = lead_shared,
            center_pos = center_pos, flank = 1000e3,
            rsid_col = "rsid", pos_col = "position",
            p_exp_col = "pval", p_out_col = "i.pval",
            exposure_label = g, outcome_label = conf$outcome)
          ggsave(file.path(dis_dir,
            sprintf("region_%s_%s_p12_%s_row%i_%s.pdf", g, gsub("\\s+","_",conf$outcome), format(p12_val), r, run_date)),
            plt, width = 8.5, height = 10)
        }
      }
    } # gene loop
  }   # p12 loop
}     # outcome loop
