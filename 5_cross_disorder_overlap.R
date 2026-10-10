# ---- Step 5 | Overlap across outcomes ---------------------------------------
# Permutation test of whether proteins were significant for several outcomes
# more often than chance (Methods):
#   - significance labels permuted within each outcome, keeping the number
#     significant fixed, 10,000 permutations
#   - the observed numbers significant in >= 2 and >= 3 outcomes compared with
#     that null
#
# Reads   mr_results_all_* from step 2
# Writes  the permutation result and the per-outcome counts behind it

library(data.table)

results_dir <- Sys.getenv("MR_OUT_DIR", "results")
run_date    <- Sys.getenv("MR_RUN_DATE", as.character(Sys.Date()))

# ---- Overlap beyond chance ----
n_perm      <- as.integer(Sys.getenv("N_PERM", "10000"))
set.seed(20260921)

files <- list.files(results_dir, paste0("^mr_results_all_.*_", run_date, "\\.csv$"), full.names = TRUE)
if (!length(files)) stop("No mr_results_all_* tables for ", run_date, " in ", results_dir)
res <- rbindlist(lapply(files, fread, select = c("id.exposure", "id.outcome", "significant")))

proteins <- sort(unique(res$id.exposure))
outcomes <- sort(unique(res$id.outcome))
count_shared <- function(sig_by_outcome) {
  n <- tabulate(match(unlist(sig_by_outcome), proteins), nbins = length(proteins))
  c(ge2 = sum(n >= 2), ge3 = sum(n >= 3))
}

tested   <- lapply(outcomes, function(o) res[id.outcome == o, id.exposure])
n_sig    <- vapply(outcomes, function(o) res[id.outcome == o, sum(significant)], integer(1))
observed <- count_shared(lapply(outcomes, function(o) res[id.outcome == o & significant == TRUE, id.exposure]))

null <- t(replicate(n_perm, count_shared(Map(function(t, k) t[sample.int(length(t), k)], tested, n_sig))))

summary_dt <- data.table(
  statistic      = c("proteins significant in >= 2 outcomes", "proteins significant in >= 3 outcomes"),
  observed       = observed,
  null_mean      = colMeans(null),
  null_95th      = apply(null, 2, quantile, 0.95),
  fold_over_null = observed / colMeans(null),
  p_empirical    = (1 + colSums(null >= rep(observed, each = n_perm))) / (1 + n_perm),
  permutations   = n_perm)
print(summary_dt)
fwrite(summary_dt, file.path(results_dir, paste0("overlap_permutation_", run_date, ".csv")))
fwrite(data.table(outcome = outcomes, proteins_tested = lengths(tested), significant = n_sig),
       file.path(results_dir, paste0("overlap_permutation_inputs_", run_date, ".csv")))
