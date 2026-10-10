# ---- Step 7 | Manuscript tables and figures ---------------------------------
# Builds Figure 2, Tables 1 and 2, the supplementary workbook (ST1-ST5) and the
# supplementary figures (SF1-SF11) from the results of steps 1c-6.
#
# Reads   step 2 results (main and the F >= 10 run), step 3a and 3b
#         colocalisation, step 6 druggability, the locus counts from 1c and
#         the UKB-PPP cis signals
# Writes  Figure 2, the tables, the workbook and SF1-SF11

library(data.table)
library(ggplot2)
library(scales)
library(ggtext)
library(patchwork)

# ---- Settings ----
env <- function(name, default) { v <- Sys.getenv(name); if (nzchar(v)) v else default }
mr_dir    <- env("MR_OUT_DIR", "results")
f10_dir   <- env("MR_F10_DIR", "results_F10")                                     # step 2 with MR_F_MIN=10
coloc_dir <- strsplit(env("COLOC_DIR", "coloc_results"), ",")[[1]]
gwas_dir  <- env("MR_GWAS_DIR", "../MR_data/target_GWASes_updated")
run_date  <- env("MR_RUN_DATE", as.character(Sys.Date()))
loci_file <- env("LOCI_FILE", file.path(gwas_dir, "outcome_loci_counts.csv"))   # from step 1c
signals_file <- env("MR_SIGNALS_FILE", "../MR_data/conditionally_independent_snps.xlsx")
out_dir   <- env("FIG_OUT_DIR", "tables_and_figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# mental illnesses first, then treatment response phenotypes; each group ordered by
# the effective sample size used in the analysis
outcomes <- c("Major depressive disorder", "Post-traumatic stress disorder", "Bipolar disorder",
              "Opioid use disorder", "Schizophrenia", "Obsessive compulsive disorder", "Panic disorder",
              "Treatment resistance in depression", "Treatment resistance in schizophrenia", "Antidepressant response")
palette_cols <- c("Major depressive disorder" = "#0072B2", "Schizophrenia" = "#D55E00",
                  "Bipolar disorder" = "#009E73", "Post-traumatic stress disorder" = "#CC79A7",
                  "Obsessive compulsive disorder" = "#E69F00", "Panic disorder" = "#56B4E9",
                  "Opioid use disorder" = "#F0E442", "Treatment resistance in depression" = "#000000",
                  "Treatment resistance in schizophrenia" = "#999999", "Antidepressant response" = "#661100")
# effect scale of each outcome's betas
scale_of <- c("Post-traumatic stress disorder" = "standardised", "Opioid use disorder" = "standardised",
              "Treatment resistance in depression" = "quantitative", "Antidepressant response" = "quantitative")
effect_scale <- function(o) ifelse(o %in% names(scale_of), scale_of[o], "log odds")

# ---- Results ----
lf <- function(pattern) unlist(lapply(coloc_dir, list.files, pattern = pattern, recursive = TRUE, full.names = TRUE))
rd <- function(d, pat) rbindlist(lapply(list.files(d, pat, full.names = TRUE), function(f) if (file.size(f) > 5) fread(f)), fill = TRUE)
mr <- rd(mr_dir, paste0("^mr_results_all_.*_", run_date, "\\.csv$"))
mr[method == "Inverse variance weighted", method := "Inverse-variance weighted"]   # hyphenated as in the text
sig <- mr[significant == TRUE, .(outcome = id.outcome, protein = id.exposure, method, nsnp, b, se, pval, p_bonferroni, Q_pval)]

# instruments in the MHC (chr6:25-34 Mb, GRCh37)
harm <- rbindlist(lapply(list.files(mr_dir, paste0("^harmonised_.*_", run_date, "\\.rds$"), full.names = TRUE), function(f) {
  x <- readRDS(f)
  rbindlist(lapply(names(x$kept), function(g) x$kept[[g]]$d[, .(
    outcome, protein = g, uniprot = protein, SNP, chr = as.integer(chr.exposure), pos = pos.exposure,
    effect_allele = effect_allele.exposure, other_allele = other_allele.exposure, eaf = eaf.exposure,
    beta.exposure, se.exposure, pval.exposure, F_stat, beta.outcome, se.outcome, pval.outcome)]))
}))
mhc <- unique(harm[chr == 6 & pos >= 25e6 & pos <= 34e6, protein])
sig[, mhc := protein %in% mhc]

# colocalisation (coloc.susie, best signal per pair) and instrument overlap
su <- rbindlist(lapply(lf("^susie_coloc_summary_.*\\.csv$"), fread), fill = TRUE)
h4 <- dcast(su[, .(h4 = max(PP.H4.abf)), by = .(outcome, protein = gene, p12)], outcome + protein ~ p12, value.var = "h4")
setnames(h4, c("5e-06", "1e-05"), c("h4_5e6", "h4_1e5"))
sig <- merge(sig, h4, by = c("outcome", "protein"), all.x = TRUE)
status <- unique(rbindlist(lapply(lf("^pair_status_.*\\.csv$"), fread)), by = c("outcome", "gene", "p12"))
abf <- dcast(status, outcome + gene ~ p12, value.var = "abf_h4")
setnames(abf, c("gene", "5e-06", "1e-05"), c("protein", "abf_5e6", "abf_1e5"))
sig <- merge(sig, abf, by = c("outcome", "protein"), all.x = TRUE)
# tier 2 is one step below tier 1: the sensitivity prior or coloc.abf, not both
sig[, tier := fcase(!is.na(h4_5e6) & h4_5e6 >= 0.8, 1L,
                    pmax(h4_1e5, abf_5e6, na.rm = TRUE) >= 0.8, 2L,
                    default = 3L)]
sig[, tier_basis := fcase(!is.na(h4_5e6) & h4_5e6 >= 0.8, "coloc.susie, main prior",
                          !is.na(h4_1e5) & h4_1e5 >= 0.8, "coloc.susie, sensitivity prior only",
                          !is.na(abf_5e6) & abf_5e6 >= 0.8, "coloc.abf, main prior",
                          !is.na(abf_1e5) & abf_1e5 >= 0.8, "coloc.abf, sensitivity prior only",
                          default = "No colocalisation")]
# instrument overlap with the credible set of the prior that gave the colocalisation
ov_files <- lf("^coloc_instrument_overlap_.*\\.csv$")
ov <- fread(ov_files[order(basename(ov_files), decreasing = TRUE)][1])
sig[, ov_p12 := fifelse(tier_basis == "coloc.susie, sensitivity prior only", 1e-5, 5e-6)]
sig <- merge(sig, ov[, .(outcome, protein = gene, ov_p12 = p12, overlap, instrument_in_cs_ids, instrument_best_ld, max_r2_instrument_cs)],
             by = c("outcome", "protein", "ov_p12"), all.x = TRUE)

# ---- Figure 2 ----
# odds ratios for case-control outcomes; exp(2 beta) approximates the OR for the
# standardised outcomes (PTSD, OUD); exp(beta) per SD for TRD (quantitative)
fig_df <- copy(sig)[, `:=`(scale = effect_scale(outcome), k = fifelse(effect_scale(outcome) == "standardised", 2, 1))]
fig_df[, `:=`(OR = exp(k * b), OR_l95 = exp(k * (b - 1.96 * se)), OR_u95 = exp(k * (b + 1.96 * se)),
              size_val = -log10(pval),
              shape_cat = fifelse(method == "Wald ratio", "Wald ratio", fifelse(Q_pval < 0.05, "IVW: *Q p* < 0.05", "IVW: *Q p* \u2265 0.05")),
              label = fifelse(mhc, paste0(protein, "<sup>\u2020</sup>"), protein))]
legend_name <- function(o) paste0(o, fifelse(effect_scale(o) == "standardised", "<br>(approximate OR)",
                                  fifelse(effect_scale(o) == "quantitative", "<br>(exp(\u03b2) per SD)", "")))
fig_df[, outcome_lab := legend_name(outcome)]
lab_levels <- legend_name(outcomes)
# panel A groups the proteins significant for more than one outcome by how many they reached
fig_df[, n_outcomes := .N, by = protein]
multi <- fig_df[n_outcomes >= 2, unique(protein)]
counts <- sort(unique(fig_df[protein %in% multi, n_outcomes]), decreasing = TRUE)
fig_df[, outcomes_grp := match(n_outcomes, counts)]
out_lab <- setNames(paste(counts, "outcomes"), seq_along(counts))
coloc_rows <- fig_df[, tier <= 2]
shown <- intersect(lab_levels, fig_df[protein %in% multi | coloc_rows, outcome_lab])
cols <- setNames(palette_cols[outcomes], lab_levels)[shown]
size_lim <- range(fig_df[protein %in% multi | coloc_rows, size_val])
shapes <- c("IVW: *Q p* < 0.05" = 17, "IVW: *Q p* \u2265 0.05" = 16, "Wald ratio" = 15)
shapes <- shapes[names(shapes) %in% fig_df[protein %in% multi | coloc_rows, shape_cat]]
shape_order <- intersect(c("Wald ratio", "IVW: *Q p* < 0.05", "IVW: *Q p* \u2265 0.05"), names(shapes))

tier_lab <- c("1" = "Tier 1", "2" = "Tier 2")
forest <- function(d, strips = NULL) {
  d <- copy(d)
  if (!is.null(strips)) d[, label := paste(label, tier, sep = "__")]
  lev <- rev(sort(unique(d$label)))
  d[, `:=`(label = factor(label, levels = lev), outcome_lab = factor(outcome_lab, levels = shown),
           shape_cat = factor(shape_cat, levels = names(shapes)))]
  # invisible points so both panels carry identical legends and patchwork can merge them
  keys <- data.table(OR = 1, label = factor(d$label[1], levels = lev), tier = d$tier[1],
                     outcome_lab = factor(rep_len(shown, max(length(shown), length(shapes))), levels = shown),
                     shape_cat = factor(rep_len(names(shapes), max(length(shown), length(shapes))), levels = names(shapes)),
                     size_val = size_lim[1])
  d[, grp := interaction(outcome_lab, shape_cat, drop = TRUE)]
  wide <- max(d$OR_u95) / min(d$OR_l95) > 8
  xb <- if (wide) c(0.25, 0.5, 1, 2, 4, 8) else c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4, 6, 8)
  xb <- xb[xb >= min(d$OR_l95) * 0.9 & xb <= max(d$OR_u95) * 1.1]
  pd <- position_dodge(width = 0.65)
  seps <- d[, .(n = uniqueN(label)), by = tier][n > 1, .(yint = seq(1.5, n - 0.5, by = 1)), by = tier]
  if (is.null(strips)) seps <- data.table(tier = NA, yint = seq_len(max(length(lev) - 1, 0)) + 0.5)
  p <- ggplot(d, aes(x = OR, y = label, colour = outcome_lab, group = grp)) +
    geom_hline(data = seps, aes(yintercept = yint), linewidth = 0.3, colour = "grey82") +
    geom_vline(xintercept = 1, linewidth = 0.6, linetype = "dashed", colour = "grey35") +
    geom_errorbarh(aes(xmin = OR_l95, xmax = OR_u95), position = pd, height = 0, linewidth = 0.55, alpha = 0.85, show.legend = FALSE) +
    geom_point(aes(size = size_val, shape = shape_cat), position = pd, alpha = 0.95) +
    geom_point(data = keys, aes(x = OR, y = label, colour = outcome_lab, shape = shape_cat, size = size_val),
               inherit.aes = FALSE, alpha = 0) +
    scale_x_log10(breaks = xb, labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.03, 0.08))) +
    scale_colour_manual(values = cols, limits = shown, drop = FALSE, name = "Psychiatric Outcome") +
    scale_shape_manual(values = shapes, limits = names(shapes), breaks = shape_order, drop = FALSE,
                       name = "Method/Heterogeneity") +
    scale_size_continuous(name = "<b>\u2212log<sub>10</sub>(<i>p</i>)</b>", range = c(1.8, 4), limits = size_lim,
                          breaks = pretty_breaks(n = 4)) +
    labs(x = "Effect per 1 SD Higher Genetically Proxied Protein Level (95% CI)", y = "Protein (Gene Symbol)") +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 3.2, shape = 16, alpha = 1)),
           shape = guide_legend(order = 2, override.aes = list(size = 3.2, colour = "black", alpha = 1)),
           size = guide_legend(order = 3, override.aes = list(shape = 16, colour = "black", alpha = 1))) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(plot.tag = element_text(face = "bold", size = 13),
          axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 8)),
          axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 8)),
          axis.text.x = element_text(size = 11, colour = "black"),
          axis.text.y = ggtext::element_markdown(size = 11, colour = "black"),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line(linewidth = 0.4, colour = "grey88"),
          legend.title = ggtext::element_markdown(size = 11, face = "bold"),
          legend.text = ggtext::element_markdown(size = 10))
  if (!is.null(strips))
    p <- p + scale_y_discrete(labels = function(x) sub("__.*", "", x)) +
      facet_grid(tier ~ ., scales = "free_y", space = "free_y", labeller = labeller(tier = strips)) +
      theme(strip.text.y = element_text(angle = 270, size = 11, face = "bold"),
            strip.background = element_rect(fill = "grey95", colour = NA))
  p
}

fig2a <- forest(copy(fig_df[protein %in% multi])[, tier := outcomes_grp], strips = out_lab)
fig2b <- forest(fig_df[coloc_rows], strips = tier_lab)
fig2a <- fig2a +
  guides(colour = guide_legend(order = 1, ncol = 1, theme = theme(legend.key.spacing.y = unit(4, "pt")),
                               override.aes = list(size = 3.2, shape = 16, alpha = 1)),
         shape = guide_legend(order = 2, ncol = 1, override.aes = list(size = 3.2, colour = "black", alpha = 1)),
         size = guide_legend(order = 3, ncol = 1, override.aes = list(shape = 16, colour = "black", alpha = 1)))
fig2 <- (fig2a | (fig2b + guides(colour = "none", shape = "none", size = "none"))) +
  plot_layout(widths = c(1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A", tag_suffix = ")") &
  theme(legend.position = "right", legend.box = "vertical", legend.box.just = "left",
        plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(out_dir, paste0("Figure2_", run_date, ".pdf")), fig2, width = 17, height = 9, device = cairo_pdf)
ggsave(file.path(out_dir, paste0("Figure2_", run_date, ".png")), fig2, width = 17, height = 9, dpi = 150)

# ---- Table 1 ----
t1 <- data.table(
  outcome = c("Major depressive disorder", "Schizophrenia", "Bipolar disorder", "Post-traumatic stress disorder",
              "Obsessive compulsive disorder", "Panic disorder", "Opioid use disorder",
              "Treatment resistance in depression", "Treatment resistance in schizophrenia", "Antidepressant response"),
  source = c("PGC MDD Working Group (2025); EUR, excluding 23andMe and UK Biobank",
             "Trubetskoy et al. (2022); European cohorts",
             "O'Connell et al. (2025); EUR, excluding UK Biobank and 23andMe",
             "Nievergelt et al. (2024); EUR",
             "Strom et al. (2025); EUR, excluding 23andMe and UK Biobank",
             "FinnGen release 12 (Kurki et al. 2023), endpoint F5_PANIC",
             "Deak et al. (2022); EUR MTAG of OUD with alcohol and cannabis use disorder",
             "Kang et al. (2024); MVP, Geisinger, Vanderbilt and Mass General Brigham",
             "Pardiñas et al. (2022); TRS vs non-TRS interaction",
             "Pain et al. (2022); percentage improvement, EUR"),
  ancestry = c("European", "European", "European", "European", "European", "Finnish", "European", "European", "European", "European"),
  n_cases = c(357636, 53386, 57833, 137136, 22717, 6911, 15251, NA, 30826, NA),
  n_controls = c(1281936, 77258, 722909, 1085746, 988884, 444414, 538935, NA, 54664, NA),
  n_total = c(1639572, 130644, 780742, 1222882, 1011601, 451325, 554186, 154506, 85490, 5218),
  # effective sample size: the largest per-variant value in each cleaned file, which is the
  # number of individuals for the two quantitative outcomes
  n_effective = c(967078, 117498, 157693, 641533, 73413, 27221, 128748, 154506, 73492, 5218),
  effect = c(rep("log odds", 3), "standardised (z-based)", "log odds", "log odds", "standardised (MTAG)",
             "per SD of predicted TRD probability", "difference in log odds (TRS vs non-TRS)", "per SD of percentage improvement"),
  covariates = c("Principal components; age and sex in some cohorts", "Principal components",
                 "Principal components (PC1-5)", "Principal components; batch and centre in UK Biobank",
                 "Principal components; age and sex in some cohorts", "Sex, age, 10 PCs, genotyping batch",
                 "Age, sex, principal components", "Age, sex, principal components (by site)",
                 "Principal components", "Age, sex, 20 principal components"),
  uk_biobank = c("No", "No", "No", "Yes", "No", "No", "No", "No", "No", "No"))
t1[outcome == "Opioid use disorder", note := "Cases and controls of the OUD GWAS; MTAG equivalent N 128,748"]
t1[outcome == "Treatment resistance in schizophrenia", note := "Schizophrenia cases: 10,501 treatment-resistant (CLOZUK) and 20,325 non-resistant (PGC 2014), each with their own controls"]
loci <- fread(loci_file)
if ("n_loci" %in% names(loci)) setnames(loci, "n_loci", "loci")
t1 <- merge(t1, loci[, .(outcome, gw_loci = loci)], by = "outcome", sort = FALSE)
t1 <- t1[match(outcomes, outcome)]
fwrite(t1, file.path(out_dir, paste0("Table1_", run_date, ".csv")))

# ---- Table 2 ----
# tiers 1 and 2; PP.H4 from coloc.susie, or from coloc.abf (marked b) where SuSiE found no credible set
t2 <- fig_df[tier <= 2]
t2[, or_ci := fifelse(scale == "quantitative", "–",
                      sprintf("%.2f (%.2f-%.2f)%s", OR, OR_l95, OR_u95, fifelse(scale == "standardised", "a", "")))]
h4_text <- function(susie, abf) fifelse(!is.na(susie), sprintf("%.2f", susie), fifelse(!is.na(abf), sprintf("%.2fb", abf), "–"))
t2[, overlap_iv := fifelse(!is.na(instrument_in_cs_ids) & instrument_in_cs_ids != "", instrument_in_cs_ids,
                   fifelse(overlap %in% TRUE, paste0(instrument_best_ld, " (r2 = ", sprintf("%.2f", max_r2_instrument_cs), ")"),
                   fifelse(is.na(overlap), "n/a", "–")))]
# results are listed by yield: within a tier the outcome with the most associations comes first,
# ties broken by the larger effective sample size
eff_n <- setNames(t1$n_effective, t1$outcome)
t2[, n_assoc := .N, by = .(tier, outcome)]
t2 <- t2[order(tier, -n_assoc, -eff_n[outcome], p_bonferroni),
         .(tier, outcome, protein, method = fifelse(method == "Wald ratio", "Wald ratio", paste0("IVW (", nsnp, ")")),
           beta_se = sprintf("%.3f (%.3f)", b, se), or_ci, p_bonferroni = signif(p_bonferroni, 2),
           Q_pval = fifelse(is.na(Q_pval), "–", as.character(signif(Q_pval, 2))),
           pph4 = h4_text(h4_5e6, abf_5e6), pph4_sensitivity = h4_text(h4_1e5, abf_1e5), overlap_iv)]
fwrite(t2, file.path(out_dir, paste0("Table2_", run_date, ".csv")))

# ---- Supplementary tables ----
f10 <- rd(f10_dir, paste0("^mr_results_all_.*_", run_date, "\\.csv$"))
f10[method == "Inverse variance weighted", method := "Inverse-variance weighted"]
sens <- rd(mr_dir, paste0("^sensitivity_.*_", run_date, "\\.csv$"))
two <- rd(mr_dir, paste0("^two_instrument_wald_ratios_.*_", run_date, "\\.csv$"))
loo <- rd(mr_dir, paste0("^leave_one_out_.*_", run_date, "\\.csv$"))
ov_all <- ov
drug <- fread(file.path(mr_dir, paste0("druggability_", run_date, ".csv")))
dgidb <- fread(file.path(mr_dir, paste0("druggability_dgidb_check_", run_date, ".csv")))

# instruments that are not UKB-PPP lead variants entered as LD proxies
signals <- as.data.table(readxl::read_excel(signals_file, skip = 4, .name_repair = "minimal"))
setnames(signals, 1, "protein_id")
signals <- signals[`Cis/trans` == "cis", .(protein = sub(":.*", "", protein_id), SNP = rsID, ukbppp_variant = `Variant ID`, PIP)]
signals <- signals[order(-PIP)][!duplicated(signals[order(-PIP), .(protein, SNP)])]
harm <- merge(harm, signals, by = c("protein", "SNP"), all.x = TRUE)
harm[, proxy := is.na(ukbppp_variant)]

sig[, sig_F10 := paste(outcome, protein) %in% f10[significant == TRUE, paste(id.outcome, id.exposure)]]

or_text <- function(b, se, o) {
  k <- fifelse(effect_scale(o) == "standardised", 2, 1)
  fifelse(effect_scale(o) == "quantitative", "",
          sprintf("%.2f (%.2f-%.2f)%s", exp(k * b), exp(k * (b - 1.96 * se)), exp(k * (b + 1.96 * se)),
                  fifelse(k == 2, "a", "")))
}
ord <- function(d) d[order(match(outcome, outcomes), p_bonferroni)]

st2 <- harm[order(match(outcome, outcomes), protein, chr, pos),
            .(`Psychiatric outcome` = outcome, `Protein (gene symbol)` = protein, UniProt = uniprot, SNP, `Chromosome` = chr,
              `Position (GRCh37)` = pos, `Effect allele` = effect_allele, `Other allele` = other_allele,
              `Effect allele frequency` = eaf, `LD proxy` = proxy,
              `UKB-PPP variant ID (GRCh38)` = ukbppp_variant, `UKB-PPP PIP` = PIP,
              `β (protein)` = beta.exposure, `SE (protein)` = se.exposure, `p (protein)` = pval.exposure,
              `F-statistic` = F_stat, `β (outcome)` = beta.outcome, `SE (outcome)` = se.outcome,
              `p (outcome)` = pval.outcome)]

iv <- harm[, .(snps = paste(SNP, collapse = "; "), mean_F = mean(F_stat), max_F = max(F_stat), n_proxy = sum(proxy)),
           by = .(outcome, protein)]
uni <- unique(harm[, .(protein, uniprot)])
full <- merge(mr[, .(outcome = id.outcome, protein = id.exposure, method, nsnp_input, n_absent_panel, n_collinear, nsnp,
                     b, se, pval, p_bonferroni, significant, Q, Q_df, Q_pval, min_F, max_r2)],
              iv, by = c("outcome", "protein"), all.x = TRUE)
full <- merge(full, uni, by = "protein", all.x = TRUE)
st3a <- merge(sig[, .(outcome, protein, mhc, h4_5e6, h4_1e5, abf_5e6, abf_1e5, tier, tier_basis, sig_F10)],
              full, by = c("outcome", "protein"))
st3a <- ord(st3a)[, .(`Psychiatric outcome` = outcome, `Protein (gene symbol)` = protein, UniProt = uniprot, `MR estimator` = method,
                      `Instruments` = nsnp, `LD proxies` = n_proxy, β = b, SE = se,
                      `OR (95% CI)` = or_text(b, se, outcome), `Raw p` = pval,
                      `Bonferroni-adjusted p` = p_bonferroni, `Cochran's Q` = Q, `Cochran's Q p-value` = Q_pval,
                      `Min. F-stat` = min_F, `Mean F-stat` = mean_F, `Max. F-stat` = max_F, SNPs = snps,
                      `MHC region` = mhc, `Significant with F ≥ 10` = sig_F10,
                      `PP.H4 (p12 = 5e-6)` = h4_5e6, `PP.H4 (p12 = 1e-5)` = h4_1e5,
                      `coloc.abf PP.H4 (p12 = 5e-6)` = abf_5e6, `coloc.abf PP.H4 (p12 = 1e-5)` = abf_1e5,
                      `Evidence tier` = tier, `Tier basis` = tier_basis)]

st3b <- sig[, .(n = .N, outcomes = paste(outcome[order(match(outcome, outcomes))], collapse = "; ")), by = protein][
  order(-n, protein), .(`Protein (gene symbol)` = protein, `n sig. outcome associations` = n, `Psychiatric outcomes` = outcomes)]

setnames(sens, c("id.exposure", "id.outcome"), c("protein", "outcome"))
st3c <- merge(sig[nsnp >= 3, .(outcome, protein, b, se, pval)], sens, by = c("outcome", "protein"))[
  order(match(outcome, outcomes), protein),
  .(`Psychiatric outcome` = outcome, `Protein (gene symbol)` = protein, Instruments = nsnp, `IVW β` = b, `IVW p` = pval,
    `MR-Egger β` = egger_b, `MR-Egger SE` = egger_se, `MR-Egger p` = egger_pval,
    `MR-Egger intercept` = egger_intercept, `Intercept SE` = egger_intercept_se, `Intercept p` = egger_intercept_pval,
    `Weighted median instruments` = wm_nsnp, `Weighted median β` = wm_b, `Weighted median SE` = wm_se,
    `Weighted median p` = wm_pval, `Leave-one-out max. p` = loo_max_pval, `Leave-one-out same sign` = loo_same_sign)]
setnames(two, c("id.exposure", "id.outcome"), c("protein", "outcome"))
st3d <- two[paste(outcome, protein) %in% sig[, paste(outcome, protein)]][order(match(outcome, outcomes), protein),
  .(`Psychiatric outcome` = outcome, `Protein (gene symbol)` = protein, SNP, `Wald ratio β` = b, SE = se, p = pval)]
setnames(loo, c("id.exposure", "id.outcome"), c("protein", "outcome"))
st3e <- loo[paste(outcome, protein) %in% sig[, paste(outcome, protein)]][order(match(outcome, outcomes), protein),
  .(`Psychiatric outcome` = outcome, `Protein (gene symbol)` = protein, `SNP left out` = SNP_left_out, β = b, SE = se, p = pval)]

best <- su[, .SD[which.max(PP.H4.abf)], by = .(outcome, protein = gene, p12)]
st4 <- merge(status[, .(outcome, protein = gene, p12, status, abf_nsnps, abf_h3, abf_h4)],
             best[, .(outcome, protein, p12, nsnps, hit1, hit2, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)],
             by = c("outcome", "protein", "p12"), all.x = TRUE)
st4 <- merge(st4, ov_all[, .(outcome, protein = gene, p12, cs95_size, instrument_in_cs, instrument_in_cs_ids,
                             max_r2_instrument_cs, instrument_best_ld, overlap)],
             by = c("outcome", "protein", "p12"), all.x = TRUE)
st4 <- merge(st4, sig[, .(outcome, protein, tier)], by = c("outcome", "protein"))
st4 <- st4[order(match(outcome, outcomes), protein, -p12),
  .(`Psychiatric outcome` = outcome, `Protein (gene symbol)` = protein, p12,
    `Result` = fcase(status == "susie_result", "coloc.susie", status == "no_credible_set", "coloc.abf (no credible set)",
                     default = "coloc.abf (coloc.susie error)"),
    `SNPs` = fifelse(status == "susie_result", nsnps, abf_nsnps),
    `Protein signal lead` = hit1, `Outcome signal lead` = hit2,
    PP.H0 = PP.H0.abf, PP.H1 = PP.H1.abf, PP.H2 = PP.H2.abf,
    PP.H3 = fifelse(status == "susie_result", PP.H3.abf, abf_h3), PP.H4 = fifelse(status == "susie_result", PP.H4.abf, abf_h4),
    `95% credible set size` = cs95_size, `Instrument in credible set` = instrument_in_cs_ids,
    `Max. r^2 instrument-credible set` = max_r2_instrument_cs, `Instrument in LD with credible set` = instrument_best_ld,
    `Instrument overlaps credible set` = overlap, `Evidence tier` = tier)]

tiers <- sig[, .(best_tier = min(tier), outcomes = paste(outcome[order(match(outcome, outcomes))], collapse = "; ")), by = protein]
drug <- merge(drug, dgidb[, .(n_dgidb = .N), by = protein], by = "protein", all.x = TRUE)
st6 <- merge(tiers, drug, by = "protein", all.x = TRUE)[order(best_tier, protein),
  .(`Protein (gene symbol)` = protein, UniProt = uniprot, `Ensembl ID` = ensembl_id, `Psychiatric outcomes` = outcomes, `Best evidence tier` = best_tier,
    `Tractability group` = tractability_group, `Small molecule` = best_sm, Antibody = best_ab, PROTAC = best_pr,
    `Other modalities` = best_oc, `Max. clinical stage` = max_clinical_stage, `Approved drugs (n)` = n_drugs_approved,
    `Clinical-stage drugs (n)` = n_drugs_clinical, `Approved drugs` = approved_drugs,
    `DGIdb counts (n)` = fifelse(is.na(n_dgidb), 0L, n_dgidb),
    `HPA location` = hpa_location)]

m10 <- f10[, .(outcome = id.outcome, protein = id.exposure, method_F10 = method, nsnp_F10 = nsnp, b_F10 = b, se_F10 = se,
               p_bonferroni_F10 = p_bonferroni, significant_F10 = significant)]
st7 <- merge(mr[, .(outcome = id.outcome, protein = id.exposure, method, nsnp, b, se, p_bonferroni, significant)],
             m10, by = c("outcome", "protein"), all = TRUE)[significant %in% TRUE | significant_F10 %in% TRUE]
st7 <- ord(st7)[, .(`Psychiatric outcome` = outcome, `Protein (gene symbol)` = protein, `MR estimator` = method, Instruments = nsnp, β = b,
                    SE = se, `Bonferroni p` = p_bonferroni, Significant = significant,
                    `MR estimator (F ≥ 10)` = method_F10, `Instruments (F ≥ 10)` = nsnp_F10, `β (F ≥ 10)` = b_F10,
                    `SE (F ≥ 10)` = se_F10, `Bonferroni p (F ≥ 10)` = p_bonferroni_F10,
                    `Significant (F ≥ 10)` = significant_F10)]

st8 <- full[order(match(outcome, outcomes), pval),
  .(`Psychiatric outcome` = outcome, `Protein (gene symbol)` = protein, UniProt = uniprot, `MR estimator` = method,
    `Instruments available` = nsnp_input, `Absent from LD panel` = n_absent_panel, `Collinear (r^2 ≥ 0.95)` = n_collinear,
    Instruments = nsnp, `LD proxies` = n_proxy, β = b, SE = se, `Raw p` = pval,
    `Bonferroni-adjusted p` = p_bonferroni, Significant = significant, `Cochran's Q` = Q, `Q df` = Q_df,
    `Cochran's Q p` = Q_pval, `Min. F-stat` = min_F, `Max. instrument r^2` = max_r2, SNPs = snps)]

captions <- c(
  ST1 = "Supplementary Table 1. Instruments used in each Mendelian randomisation test, after LD-proxy substitution and Steiger filtering. Protein effects are per standard deviation of protein level, and outcome effects are on the scale of each outcome GWAS (Table 1). F = (β/SE)^2. 'LD proxy' is 1 where a proxy replaced an instrument missing from the outcome GWAS. PIP is the posterior inclusion probability from UK Biobank Pharma Proteomics Project (UKB-PPP) fine-mapping, given for instruments that are UKB-PPP lead variants rather than LD proxies.",
  ST2A = "Supplementary Table 2A. Mendelian randomisation (MR) results for all 16,621 protein-outcome tests. An instrument was excluded from a protein's estimate if it was absent from the LD reference panel, or collinear (r^2 ≥ 0.95) with a stronger instrument.",
  ST2B = "Supplementary Table 2B. Mendelian randomisation (MR) results for protein-outcome associations significant after Bonferroni correction. Odds ratios for post-traumatic stress disorder and opioid use disorder are approximate (exp(2β), marked a), and none are given for treatment resistance in depression. 'Significant with F ≥ 10' indicates whether the association remains significant after excluding instruments with F < 10. PP.H4 is the posterior probability of a shared causal variant from coloc.susie or, where it found no credible set, coloc.abf, at p12 = 5e-6 (main) and 1e-5 (sensitivity). Evidence tiers are defined in the colocalisation section of the Methods. MHC: major histocompatibility complex.",
  ST2C = "Supplementary Table 2C. Psychiatric outcomes associated with each Mendelian randomisation-significant protein.",
  ST3A = "Supplementary Table 3A. Sensitivity analyses for Mendelian randomisation-significant associations with ≥ 3 instruments. Inverse-variance weighted (IVW) estimates are from the main analysis. The weighted median was estimated after pruning instruments to r^2 < 0.10, and is given where at least three remained.",
  ST3B = "Supplementary Table 3B. Leave-one-out estimates for Mendelian randomisation-significant associations with ≥ 3 instruments.",
  ST3C = "Supplementary Table 3C. Wald ratio of each instrument, for Mendelian randomisation-significant associations with two instruments.",
  ST3D = "Supplementary Table 3D. Associations significant in the main Mendelian randomisation (MR) analysis, in the analysis excluding instruments with F < 10, or in both, with results from each analysis. Empty cells mean the protein had no instrument with F ≥ 10 and could not be tested in that analysis.",
  ST4 = "Supplementary Table 4. Colocalisation of all Mendelian randomisation-significant associations at p12 = 5e-6 (main) and 1e-5 (sensitivity). coloc.abf, which assumes a single causal variant, was used where coloc.susie did not return a result. An instrument overlaps the colocalised signal if it lies in the signal's 95% credible set or has r^2 ≥ 0.80 with a variant in it.",
  ST5 = "Supplementary Table 5. Druggability of Mendelian randomisation-significant proteins, annotated from Open Targets Platform release 26.06 (tractability and drugs) and the Human Protein Atlas (HPA; protein location). Tractability grouping follows Dardani et al. (2025), and only drugs with a recorded mechanism of action on the protein are listed. Protein location is the HPA secretome classification for secreted proteins and its predicted protein class otherwise. The DGIdb column counts approved-drug interactions listed by DGIdb 5.0 but not by Open Targets; none has such a mechanism. PROTAC: proteolysis-targeting chimera.")
# numbered in the order the manuscript first cites them: MR results (2), sensitivity analyses (3),
# colocalisation (4), druggability (5)
sheets <- list(ST1 = st2, ST2A = st8, ST2B = st3a, ST2C = st3b, ST3A = st3c, ST3B = st3e,
               ST3C = st3d, ST3D = st7, ST4 = st4, ST5 = st6)
caption_width <- 150
header_width <- 18
pad <- 0.71

wrap_lines <- function(s, chars) {
  n <- 1L; cur <- 0L
  for (word in strsplit(s, " ")[[1]]) {
    k <- nchar(word)
    if (cur == 0L) cur <- k
    else if (cur + 1L + k <= chars) cur <- cur + 1L + k
    else { n <- n + 1L; cur <- k }
  }
  n
}

wb <- openxlsx::createWorkbook()
for (nm in names(sheets)) {
  d <- sheets[[nm]]
  need <- nchar(names(d)) * 1.08 + 1
  w <- c(pmin(pmax(need, 8), header_width), rep(8.43, 40))
  cum <- cumsum(w)
  span <- if (cum[1] + pad >= caption_width) 1L else max(which(cum + pad * seq_along(cum) <= caption_width))
  w[seq_len(span)] <- w[seq_len(span)] +
    (caption_width - pad * span - cum[span]) / span
  w <- w[seq_len(max(ncol(d), span))]
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, captions[[nm]], startRow = 1)
  openxlsx::writeData(wb, nm, d, startRow = 3,
                      headerStyle = openxlsx::createStyle(textDecoration = "bold", wrapText = TRUE, valign = "bottom"))
  openxlsx::setColWidths(wb, nm, cols = seq_along(w), widths = w)
  openxlsx::mergeCells(wb, nm, cols = 1:span, rows = 1)
  openxlsx::addStyle(wb, nm, openxlsx::createStyle(wrapText = TRUE, halign = "left", valign = "center"),
                     rows = 1, cols = 1)
  openxlsx::setRowHeights(wb, nm, rows = c(1, 3), heights = c(
    15 * wrap_lines(captions[[nm]], caption_width * 1.4) + 6,
    15 * max(vapply(names(d), wrap_lines, integer(1), chars = header_width)) + 6))
}
xlsx_file <- file.path(out_dir, paste0("Supplementary_Tables_", run_date, ".xlsx"))
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)

rich_runs <- function(s, bold_chars) {
  ch <- strsplit(s, "")[[1]]; n <- length(ch)
  ital <- logical(n)
  m <- gregexpr("(?<![A-Za-z0-9])[pFQnr](?![A-Za-z])", s, perl = TRUE)[[1]]
  if (m[1] > 0) ital[m] <- TRUE
  for (fn in c("coloc.susie", "coloc.abf")) {
    at <- gregexpr(fn, s, fixed = TRUE)[[1]]
    if (at[1] > 0) for (i in at) ital[i:(i + nchar(fn) - 1L)] <- TRUE
  }
  bold <- seq_len(n) <= bold_chars
  key <- paste(bold, ital)
  grp <- cumsum(c(TRUE, key[-1] != key[-n]))
  esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x)))
  runs <- vapply(split(seq_len(n), grp), function(ix) {
    pr <- paste0(if (bold[ix[1]]) "<b/>" else "", if (ital[ix[1]]) "<i/>" else "")
    sprintf("<r>%s<t xml:space=\"preserve\">%s</t></r>", if (nzchar(pr)) paste0("<rPr>", pr, "</rPr>") else "",
            esc(paste(ch[ix], collapse = "")))
  }, character(1))
  paste0("<si>", paste(runs, collapse = ""), "</si>")
}
body_font <- "Aptos Narrow"
format_strings <- function(path, bold_chars) {
  dir <- file.path(tempdir(), "xlsx_rich"); unlink(dir, recursive = TRUE); dir.create(dir)
  utils::unzip(path, exdir = dir)
  for (part in c(file.path(dir, "xl", "styles.xml"), file.path(dir, "xl", "theme", "theme1.xml"))) {
    y <- paste(readLines(part, warn = FALSE), collapse = "")
    y <- gsub('val="Calibri"', sprintf('val="%s"', body_font), y, fixed = TRUE)
    y <- gsub('typeface="Calibri"', sprintf('typeface="%s"', body_font), y, fixed = TRUE)
    writeLines(y, part)
  }
  f <- file.path(dir, "xl", "sharedStrings.xml")
  x <- paste(readLines(f, warn = FALSE), collapse = "")
  si <- regmatches(x, gregexpr("<si>.*?</si>", x))[[1]]
  regmatches(x, gregexpr("<si>.*?</si>", x)) <- list(vapply(si, function(b) {
    t <- sub("^<si><t[^>]*>(.*)</t></si>$", "\\1", b)
    if (identical(t, b)) return(b)
    s <- gsub("&amp;", "&", gsub("&gt;", ">", gsub("&lt;", "<", gsub("&apos;", "'", t))))
    if (!s %in% names(bold_chars)) return(b)
    rich_runs(s, bold_chars[[s]])
  }, character(1)))
  writeLines(x, f)
  zip::zip(normalizePath(path), list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE),
           root = dir, mode = "mirror")
}
title_end <- function(s) { m <- gregexpr("\\. ", s)[[1]]; if (length(m) < 2 || m[1] < 0) nchar(s) else m[2] }
hdrs <- unique(unlist(lapply(sheets, names), use.names = FALSE))
fn_cells <- unique(unlist(lapply(sheets, function(d)
  unlist(lapply(d, function(col) if (is.character(col)) grep("coloc[.]", col, value = TRUE))))))
format_strings(xlsx_file, c(setNames(vapply(captions, title_end, numeric(1)), captions),
                            setNames(rep(Inf, length(hdrs)), hdrs),
                            setNames(rep(0, length(fn_cells)), fn_cells)))

# ---- Supplementary figures ----
# SF1-SF10: Manhattan plots, one per outcome, Wald ratio and IVW pooled as in the Bonferroni correction
pos <- harm[, .(chr = chr[1], pos = min(pos)), by = .(outcome, protein)]
man <- merge(mr[, .(outcome = id.outcome, protein = id.exposure, method, pval, significant, n_tests)], pos,
             by = c("outcome", "protein"))
chr_len <- man[, .(len = max(pos)), by = chr][order(chr)][, start := cumsum(as.numeric(len)) - len]
man <- merge(man, chr_len[, .(chr, start)], by = "chr")[, gpos := pos + start]
axis_df <- man[, .(centre = (min(gpos) + max(gpos)) / 2), by = chr][order(chr)]

man_shapes <- c("Wald ratio" = 16, "Inverse-variance weighted" = 17)
for (i in seq_along(outcomes)) {
  o <- outcomes[i]; d <- man[outcome == o]
  d[, colour_key := fifelse(significant, "sig", fifelse(chr %% 2 == 1, "odd", "even"))]
  thr <- -log10(0.05 / d$n_tests[1])
  p <- ggplot(d, aes(gpos, -log10(pval))) +
    geom_point(aes(colour = colour_key, shape = method), size = 1.4, alpha = 0.85) +
    geom_hline(yintercept = thr, linetype = "dashed", linewidth = 0.4) +
    ggrepel::geom_label_repel(data = d[significant == TRUE], aes(label = protein), size = 3, max.overlaps = Inf,
                              label.size = 0.15, min.segment.length = 0, box.padding = 0.3, seed = 1) +
    scale_colour_manual(values = c(odd = "grey65", even = "#9ECAE1", sig = palette_cols[[o]]), guide = "none") +
    scale_shape_manual(values = man_shapes, breaks = names(man_shapes), name = NULL) +
    scale_x_continuous(breaks = axis_df$centre, labels = axis_df$chr, expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = o, x = "Chromosome", y = expression(-log[10](italic(p)))) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
          axis.text.x = element_text(size = 8))
  ggsave(file.path(out_dir, sprintf("SF%d_Manhattan_%s_%s.pdf", i, gsub("[^A-Za-z]+", "_", o), run_date)), p,
         width = 10, height = 5, device = cairo_pdf)
}

# SF11: scatter plots of every significant multi-instrument association
multi_iv <- ord(sig[nsnp >= 2])
sc <- merge(harm, multi_iv[, .(outcome, protein, b, rank = .I)], by = c("outcome", "protein"))
sc[beta.exposure < 0, `:=`(beta.exposure = -beta.exposure, beta.outcome = -beta.outcome)]   # orient to the protein-increasing allele
sc[, panel := factor(paste0(protein, "\n", outcome), levels = multi_iv[, paste0(protein, "\n", outcome)])]
sf11 <- ggplot(sc, aes(beta.exposure, beta.outcome)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_errorbar(aes(ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), width = 0, colour = "grey50") +
  geom_errorbarh(aes(xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), height = 0, colour = "grey50") +
  geom_point(aes(shape = proxy), size = 1.5) +
  geom_abline(aes(intercept = 0, slope = b), colour = "#D55E00", linewidth = 0.6) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("UKB-PPP lead variant", "LD proxy"), name = NULL) +
  scale_x_continuous(n.breaks = 4) + scale_y_continuous(n.breaks = 4) +
  facet_wrap(~ panel, scales = "free", ncol = 5) +
  labs(x = "SNP effect on protein (SD)", y = "SNP effect on outcome") +
  theme_bw(base_size = 8.5) + theme(legend.position = "bottom", strip.text = element_text(size = 6.8))
ggsave(file.path(out_dir, paste0("SF11_scatter_plots_", run_date, ".pdf")), sf11, width = 10.5, height = 1.6 * ceiling(nrow(multi_iv) / 5) + 0.6,
       device = cairo_pdf)
