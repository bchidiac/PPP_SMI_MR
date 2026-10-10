# ---- Step 3b | Instruments in the colocalised signal ------------------------
# For each colocalising signal (PP.H4 >= 0.80, from coloc.susie or from
# coloc.abf where SuSiE found no credible set), is an MR instrument for that
# pair in the signal's 95% credible set, or in LD with it (r2 >= 0.80, EUR
# panel)? (Supplementary Methods SM5.)
#
# Reads   step 3a summaries and per-SNP posteriors, the harmonised instruments
#         from step 2, 1b LD panel
# Writes  one row per colocalising signal
# Feeds   7 (Table 2, Supplementary Table 4)

suppressMessages(library(data.table))
lf <- function(pattern) unlist(lapply(coloc_dir, list.files, pattern = pattern, recursive = TRUE, full.names = TRUE))

# ---- Settings ----
env <- function(name, default) { v <- Sys.getenv(name); if (nzchar(v)) v else default }
coloc_dir <- strsplit(env("COLOC_DIR", "coloc_results"), ",")[[1]]
mr_dir    <- env("MR_OUT_DIR", "results")
mr_date   <- env("MR_RUN_DATE", as.character(Sys.Date()))
ld_panel  <- env("MR_LD_PANEL", "../MR_data/LD_reference_panel_1000G_EUR_phase3/EUR")
out_file  <- env("OVERLAP_OUT", file.path(coloc_dir[1], paste0("coloc_instrument_overlap_", mr_date, ".csv")))
thr_h4    <- 0.80
r2_min    <- 0.80

# genotypes from the PLINK panel, for LD
bim <- fread(paste0(ld_panel, ".bim"), header = FALSE, select = c(2), col.names = "SNP")[, idx := .I]
n_ind <- nrow(fread(paste0(ld_panel, ".fam"), header = FALSE))
read_geno <- function(snps) {
  b <- bim[SNP %in% snps][!duplicated(SNP)]
  if (!nrow(b)) return(NULL)
  nb <- ceiling(n_ind / 4); con <- file(paste0(ld_panel, ".bed"), "rb"); on.exit(close(con))
  G <- t(vapply(b$idx, function(i) {
    seek(con, 3 + (i - 1) * nb); x <- as.integer(readBin(con, "raw", nb))
    c(2, NA, 1, 0)[as.vector(rbind(x %% 4, (x %/% 4) %% 4, (x %/% 16) %% 4, (x %/% 64) %% 4))[1:n_ind] + 1]
  }, numeric(n_ind)))
  rownames(G) <- b$SNP; G
}

# instruments entering the final MR test for each pair
harm <- list.files(mr_dir, paste0("^harmonised_.*_", mr_date, "\\.rds$"), full.names = TRUE)
instr <- rbindlist(lapply(harm, function(f) {
  x <- readRDS(f)
  rbindlist(lapply(names(x$kept), function(g) {
    d <- x$kept[[g]]$d
    data.table(outcome = d$outcome[1], gene = g, SNP = d$SNP)
  }))
}))

# ---- Overlap ----
summ <- rbindlist(lapply(lf("^susie_coloc_summary_.*\\.csv$"),
                         fread), fill = TRUE)
summ[, row := seq_len(.N), by = .(outcome, gene, p12)]
sig <- summ[PP.H4.abf >= thr_h4]

overlap <- function(o, g, p12, signal, h4, f) {
  if (!length(f)) return(data.table(outcome = o, gene = g, p12 = p12, signal = signal,
                                    PP.H4 = h4, note = "per-SNP posteriors not found"))
  pp <- fread(f[1])[order(-pp)]
  pp[, share := pp / sum(pp)]
  cs <- pp[seq_len(which(cumsum(share) >= 0.95)[1]), snp]            # smallest set reaching 95%
  iv <- instr[outcome == o & gene == g, SNP]
  in_cs <- intersect(iv, cs)
  max_r2 <- NA_real_; best_iv <- NA_character_
  G <- read_geno(unique(c(iv, cs)))
  if (!is.null(G)) {
    ivg <- intersect(iv, rownames(G)); csg <- intersect(cs, rownames(G))
    if (length(ivg) && length(csg)) {
      r2 <- suppressWarnings(cor(t(G[ivg, , drop = FALSE]), t(G[csg, , drop = FALSE]), use = "pairwise.complete.obs"))^2
      max_r2 <- max(r2, na.rm = TRUE)
      best_iv <- ivg[which(r2 == max_r2, arr.ind = TRUE)[1, 1]]
    }
  }
  data.table(outcome = o, gene = g, p12 = p12, signal = signal, PP.H4 = h4,
             cs95_size = length(cs), lead_shared_variant = cs[1], n_instruments = length(iv),
             instrument_in_cs = length(in_cs) > 0,
             instrument_in_cs_ids = paste(in_cs, collapse = ";"),
             max_r2_instrument_cs = max_r2, instrument_best_ld = best_iv,
             overlap = length(in_cs) > 0 | isTRUE(max_r2 >= r2_min),
             note = "")
}
out_susie <- rbindlist(lapply(seq_len(nrow(sig)), function(k) {
  s <- sig[k]
  f <- lf(sprintf("^susie_coloc_snps_%s_%s_p12_%s_row%d_.*\\.csv$", s$gene,
                  gsub("\\s+", "_", s$outcome), format(s$p12), s$row))
  overlap(s$outcome, s$gene, s$p12, as.character(s$row), s$PP.H4.abf, f)
}), fill = TRUE)

# coloc.abf fallback pairs (no SuSiE credible set, or SuSiE did not converge)
status <- rbindlist(lapply(lf("^pair_status_task.*\\.csv$"), fread))
abf <- unique(status[status != "susie_result" & abf_h4 >= thr_h4], by = c("outcome", "gene", "p12"))
out_abf <- rbindlist(lapply(seq_len(nrow(abf)), function(k) {
  s <- abf[k]
  f <- lf(sprintf("^abf_coloc_snps_%s_%s_p12_%s_.*\\.csv$", s$gene,
                  gsub("\\s+", "_", s$outcome), format(s$p12)))
  overlap(s$outcome, s$gene, s$p12, "abf", s$abf_h4, f)
}), fill = TRUE)
out <- rbind(out_susie, out_abf, fill = TRUE)

fwrite(out, out_file)
for (p in sort(unique(out$p12))) {
  o <- out[p12 == p]
  cat(sprintf("p12 %g: colocalising signals %d (pairs %d); instrument in or r2 >= %.2f with the 95%% credible set: %d signals\n",
              p, nrow(o), uniqueN(o[, paste(outcome, gene)]), r2_min, sum(o$overlap, na.rm = TRUE)))
}
