## ---------------------------------------------------------------------------
## 01_run_pilot.R -- small Monte Carlo (R = 5 per scenario, smaller N) used
## to validate the simulation plumbing and to spot bugs before the full run.
##
## Output: results/pilot_metrics.rds  (long data frame of per-replicate
## metrics for each method x scenario)
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

PILOT <- modifyList(DEFAULT_DESIGN, list(
  N = 400, V = 150, K = 5, P = 4,
  R = 5,
  boot_B = 80,           ## small bootstrap for pilot
  fit_max_iter = 50,
  boot_max_iter = 20
))

set.seed(2026)
master_seeds <- sample.int(.Machine$integer.max, PILOT$R)

all_rows <- list()
t_pilot <- Sys.time()

for (sc in PILOT$scenarios) {
  cat("\n==== scenario:", sc, "====\n")
  for (r in seq_len(PILOT$R)) {
    cat(sprintf("  replicate %d/%d  (seed %d) ...", r, PILOT$R, master_seeds[r]))
    t0 <- Sys.time()
    sim <- simulate_gscamm(N = PILOT$N, V = PILOT$V, K = PILOT$K, P = PILOT$P,
                           scenario = sc, seed = master_seeds[r],
                           doc_length_mean = if (sc == "high_sparsity") 60 else 300)
    rows <- tryCatch(run_one_replicate(sim, sc, design = PILOT),
                     error = function(e) {
                       message(" FAILED: ", conditionMessage(e)); NULL
                     })
    if (!is.null(rows)) {
      rows$replicate <- r
      rows$seed <- master_seeds[r]
      all_rows[[length(all_rows) + 1L]] <- rows
    }
    cat(sprintf(" %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))
  }
}

results <- do.call(rbind, all_rows)
results$method <- factor(results$method,
                         levels = c("gscamm", "gscamm_boot", "lda", "stm"))

cat("\n==== pilot summary ====\n")
agg <- aggregate(
  cbind(rmse_theta, rmse_phi, rmse_B, coverage_B, width_B, perplexity, time)
    ~ scenario + method,
  data = results,
  FUN = function(x) round(mean(x, na.rm = TRUE), 3)
)
print(agg, row.names = FALSE)

saveRDS(results, file.path(RESULT_DIR, "pilot_metrics.rds"))
write.csv(results, file.path(RESULT_DIR, "pilot_metrics.csv"), row.names = FALSE)
cat(sprintf("\nPilot finished in %.1f s. Results saved to %s\n",
            as.numeric(Sys.time() - t_pilot, units = "secs"),
            file.path(RESULT_DIR, "pilot_metrics.rds")))
