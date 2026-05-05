## ---------------------------------------------------------------------------
## 01b_run_validation.R -- validation pilot for the parametric bootstrap.
##
## Confirms that under the deterministic GSCA-MM model the parametric
## bootstrap (W^(b) ~ Mult(L, theta_hat' Phi_hat) + warm-started refit)
## delivers near-nominal CI coverage on B, in contrast to:
##   - the plug-in WLS sandwich (under-covers, ~0.70)
##   - the legacy noise-augmented row bootstrap (over-covers, ~0.99)
##
## Setup: paper-scale design (N=1000, V=500, K=10, P=8) at R=30 replicates
## per scenario, B=200 bootstrap reps, all 3 scenarios, the three model
## comparators (gscamm, lda, stm). Stops before the production R=100 run
## so the numbers can be reviewed first.
##
## Outputs:
##   results/validation_metrics.csv   per-replicate metrics (long format)
##   results/validation_table.csv     aggregated mean/sd table
##   stdout                          a summary of the three coverages
##                                   (plug-in / parametric / legacy noise)
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

LOG_FILE <- file.path(RESULT_DIR, "validation_pilot.log")
.log_con <- file(LOG_FILE, open = "wt")
.log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(...), "\n")
  cat(msg); cat(msg, file = .log_con); flush(.log_con)
}

R_arg <- as.integer(Sys.getenv("GSCAMM_VALIDATION_R", unset = "30"))
B_arg <- as.integer(Sys.getenv("GSCAMM_VALIDATION_B", unset = "200"))

DESIGN <- modifyList(DEFAULT_DESIGN, list(
  N = 1000, V = 500, K = 10, P = 8,
  R = R_arg,
  boot_B = B_arg,                ## legacy noise-augmented bootstrap
  param_boot_B = B_arg,          ## new parametric bootstrap (the focus)
  fit_max_iter = 200,
  boot_max_iter = 60,
  param_boot_max_iter = 50,
  noise_scale = 1.0,
  use_polish = TRUE,
  sigma2_polish = 0.25,
  use_stm_spectral = FALSE
))

.log("=== Validation pilot for parametric bootstrap ===")
.log(sprintf("  R = %d replicates per scenario", DESIGN$R))
.log(sprintf("  B = %d (parametric)  B = %d (legacy noise-aug)",
             DESIGN$param_boot_B, DESIGN$boot_B))
.log(sprintf("  scenarios = %s",
             paste(DESIGN$scenarios, collapse = ", ")))

n_phys <- max(1L, parallel::detectCores(logical = FALSE))
n_cores <- max(1L, n_phys - 1L)
if (.Platform$OS.type == "windows") n_cores <- 1L
.log(sprintf("  cores = %d (OS: %s)", n_cores, .Platform$OS.type))

doc_lengths <- list(baseline = 1000, high_covariate = 1000,
                    high_sparsity = 20)

set.seed(2026)
master_seeds <- sample.int(.Machine$integer.max, DESIGN$R)

.run_worker <- function(r, sc, design, doc_lengths, master_seeds) {
  t0 <- Sys.time()
  sim <- simulate_gscamm(N = design$N, V = design$V,
                         K = design$K, P = design$P,
                         scenario = sc, seed = master_seeds[r],
                         doc_length_mean = doc_lengths[[sc]])
  rows <- tryCatch(run_one_replicate(sim, sc, design = design),
                   error = function(e) {
                     message(sprintf("[%s rep %d] FAILED: %s",
                                     sc, r, conditionMessage(e)))
                     NULL
                   })
  if (!is.null(rows)) {
    rows$replicate <- r
    rows$seed <- master_seeds[r]
  }
  list(rows = rows, r = r,
       elapsed = as.numeric(Sys.time() - t0, units = "secs"))
}

chunk_ids <- function(R, size) split(seq_len(R), ceiling(seq_len(R) / size))

all_rows   <- list()
t_global   <- Sys.time()
total_done <- 0L
total_reps <- length(DESIGN$scenarios) * DESIGN$R

for (sc in DESIGN$scenarios) {
  t_sc <- Sys.time()
  .log(sprintf("========== scenario: %s [R=%d, %d cores] ==========",
               sc, DESIGN$R, n_cores))
  sc_rows <- list()
  for (chunk in chunk_ids(DESIGN$R, n_cores)) {
    res_chunk <- if (n_cores > 1L)
      parallel::mclapply(chunk, .run_worker,
                         sc = sc, design = DESIGN,
                         doc_lengths = doc_lengths,
                         master_seeds = master_seeds,
                         mc.cores = n_cores, mc.preschedule = FALSE)
    else
      lapply(chunk, .run_worker,
             sc = sc, design = DESIGN,
             doc_lengths = doc_lengths,
             master_seeds = master_seeds)
    for (res in res_chunk) {
      total_done <- total_done + 1L
      pct <- 100 * total_done / total_reps
      elapsed_g <- as.numeric(Sys.time() - t_global, units = "mins")
      .log(sprintf("  %s rep %3d/%d  %6.1fs  | total: %5.1f min  [%5.1f%%]",
                   sc, res$r, DESIGN$R, res$elapsed, elapsed_g, pct))
      if (!is.null(res$rows))
        sc_rows[[length(sc_rows) + 1L]] <- res$rows
    }
  }
  all_rows <- c(all_rows, sc_rows)
  .log(sprintf("  --> scenario %s complete in %.1f min",
               sc, as.numeric(Sys.time() - t_sc, units = "mins")))
}

results <- do.call(rbind, all_rows)
results$method <- factor(
  results$method,
  levels = intersect(c("gscamm", "gscamm_boot", "lda", "stm", "stm_spectral"),
                     unique(as.character(results$method)))
)

saveRDS(results, file.path(RESULT_DIR, "validation_metrics.rds"))
write.csv(results, file.path(RESULT_DIR, "validation_metrics.csv"),
          row.names = FALSE)

## ----- aggregated table -------------------------------------------------
.agg <- function(x) round(mean(x, na.rm = TRUE), 4)
agg <- aggregate(
  cbind(rmse_theta, rmse_theta_map, rmse_B,
        coverage_B, width_B,
        coverage_B_boot_param, width_B_boot_param,
        time, time_param_boot)
    ~ scenario + method, data = results, FUN = .agg, na.action = na.pass)
write.csv(agg, file.path(RESULT_DIR, "validation_table.csv"),
          row.names = FALSE)

## ----- terminal summary -------------------------------------------------
.log("\n=== Coverage and width summary ===")
.log("Each cell shows the mean over R replicates (B=", DESIGN$param_boot_B, ").")
for (sc in DESIGN$scenarios) {
  .log(sprintf("\n--- %s ---", sc))
  pi_row <- agg[agg$scenario == sc & agg$method == "gscamm", ]
  bo_row <- agg[agg$scenario == sc & agg$method == "gscamm_boot", ]
  if (nrow(pi_row)) {
    .log(sprintf("  plug-in WLS         coverage=%.3f  width=%.3f  (target 0.95)",
                 pi_row$coverage_B, pi_row$width_B))
    .log(sprintf("  parametric boot     coverage=%.3f  width=%.3f  (target 0.93-0.97)",
                 pi_row$coverage_B_boot_param,
                 pi_row$width_B_boot_param))
  }
  if (nrow(bo_row)) {
    .log(sprintf("  legacy noise-aug    coverage=%.3f  width=%.3f  (was ~0.99)",
                 bo_row$coverage_B, bo_row$width_B))
  }
}

total_min <- as.numeric(Sys.time() - t_global, units = "mins")
.log(sprintf("\n=== Validation pilot complete in %.1f min ===", total_min))
.log(paste("metrics ->",
           file.path(RESULT_DIR, "validation_metrics.{rds,csv}")))
.log(paste("table   ->", file.path(RESULT_DIR, "validation_table.csv")))
.log(paste("log     ->", LOG_FILE))
.log("\nNext: review the parametric coverage; if in [0.90, 0.99] proceed to")
.log("the production run via Rscript code/02_run_full.R.")
close(.log_con)
