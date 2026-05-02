## ---------------------------------------------------------------------------
## 02_run_full.R -- full Monte Carlo, paper setup
##   N = 1000, V = 500, K = 10, P = 8, R = 100 per scenario
##
## Outputs:
##   results/full_metrics.rds  (long data frame)
##   results/full_metrics.csv
## ---------------------------------------------------------------------------

.find_setup <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).*", args, perl = TRUE))
  here <- if (length(m)) dirname(normalizePath(m[1])) else getwd()
  if (!file.exists(file.path(here, "00_setup.R")))
    here <- file.path(getwd(), "replication_package/simulations/code")
  file.path(here, "00_setup.R")
}
source(.find_setup())

## allow override via environment variable for ad-hoc smaller runs
R_arg <- Sys.getenv("GSCAMM_SIM_R", unset = "100")
B_arg <- Sys.getenv("GSCAMM_BOOT_B", unset = "200")

DESIGN <- modifyList(DEFAULT_DESIGN, list(
  N = 1000, V = 500, K = 10, P = 8,
  R = as.integer(R_arg),
  boot_B = as.integer(B_arg)
))
cat("=== Full simulation design ===\n")
str(DESIGN)

## per-scenario doc-length, matching the paper
doc_lengths <- list(baseline       = 1000,
                    high_covariate = 1000,
                    high_sparsity  = 20)

set.seed(2026)
master_seeds <- sample.int(.Machine$integer.max, DESIGN$R)

all_rows <- list()
t_global <- Sys.time()

for (sc in DESIGN$scenarios) {
  cat("\n=========== scenario:", sc, "===========\n")
  for (r in seq_len(DESIGN$R)) {
    t0 <- Sys.time()
    sim <- simulate_gscamm(N = DESIGN$N, V = DESIGN$V, K = DESIGN$K, P = DESIGN$P,
                           scenario = sc, seed = master_seeds[r],
                           doc_length_mean = doc_lengths[[sc]])
    rows <- tryCatch(run_one_replicate(sim, sc, design = DESIGN),
                     error = function(e) {
                       message(sprintf("[%s rep %d] FAILED: %s",
                                       sc, r, conditionMessage(e)))
                       NULL
                     })
    if (!is.null(rows)) {
      rows$replicate <- r
      rows$seed <- master_seeds[r]
      all_rows[[length(all_rows) + 1L]] <- rows
    }
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    total_elapsed <- as.numeric(Sys.time() - t_global, units = "mins")
    cat(sprintf("  %s rep %3d/%d done in %5.1fs   total: %5.1f min\n",
                sc, r, DESIGN$R, elapsed, total_elapsed))
    ## checkpoint every 10 replicates per scenario
    if (r %% 10L == 0L) {
      results <- do.call(rbind, all_rows)
      saveRDS(results, file.path(RESULT_DIR, "full_metrics_partial.rds"))
    }
  }
}

results <- do.call(rbind, all_rows)
results$method <- factor(results$method,
                         levels = c("gscamm", "gscamm_boot", "lda", "stm"))

saveRDS(results, file.path(RESULT_DIR, "full_metrics.rds"))
write.csv(results, file.path(RESULT_DIR, "full_metrics.csv"), row.names = FALSE)
unlink(file.path(RESULT_DIR, "full_metrics_partial.rds"))

cat(sprintf("\n=== Full run complete in %.1f min ===\n",
            as.numeric(Sys.time() - t_global, units = "mins")))
cat("Results saved to:\n  ", file.path(RESULT_DIR, "full_metrics.rds"), "\n")
