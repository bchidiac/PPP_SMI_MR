# ---- Step 4 | Estimator use and heterogeneity -------------------------------
# How often each MR estimator was used across outcomes, and the share of IVW
# tests with between-instrument heterogeneity (Cochran's Q p < 0.05).
#
# Reads   plot_raw_* and heterogeneity_stats_* from step 2
# Writes  two summary tables
# Option  MR_RUN_DATE selects which run's tables to summarise

library(data.table)
library(dplyr)

results_dir <- Sys.getenv("MR_OUT_DIR", "results")
run_date    <- Sys.getenv("MR_RUN_DATE", as.character(Sys.Date()))

# ---- Estimator usage ----
raw_files <- list.files(results_dir,
                        pattern = paste0("^plot_raw_.*_", run_date, "\\.csv$"),
                        full.names = TRUE)

raw_dt <- rbindlist(lapply(raw_files, fread), fill = TRUE) %>%
  mutate(outcome = ifelse(!is.na(outcome), outcome, id.outcome)) %>%
  select(outcome, protein = id.exposure, method)

method_summary <- raw_dt %>%
  distinct(outcome, protein, method) %>%
  count(method, name = "n_tests") %>%
  mutate(prop = n_tests / sum(n_tests))

print(method_summary)
fwrite(as.data.table(method_summary),
       file.path(results_dir, paste0("estimator_usage_", run_date, ".csv")))

# ---- Heterogeneity ----
het_files <- list.files(results_dir,
                        pattern = paste0("^heterogeneity_stats_.*_", run_date, "\\.csv$"),
                        full.names = TRUE)

if (length(het_files) > 0) {
  het_dt <- rbindlist(lapply(het_files, fread), fill = TRUE) %>%
    filter(grepl("Inverse variance weighted|Inverse_variance_weighted|IVW", method)) %>%
    rename(protein = id.exposure) %>%
    mutate(outcome = ifelse(!is.na(outcome), outcome, id.outcome))

  het_summary <- het_dt %>%
    distinct(outcome, protein, Q_pval) %>%
    summarise(
      n_ivw_tests_with_Q = n(),
      n_Qp_lt_0.05       = sum(Q_pval < 0.05, na.rm = TRUE),
      prop_Qp_lt_0.05    = n_Qp_lt_0.05 / n_ivw_tests_with_Q)

  print(het_summary)
  fwrite(as.data.table(het_summary),
         file.path(results_dir, paste0("heterogeneity_summary_", run_date, ".csv")))
}
