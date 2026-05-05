## ---------------------------------------------------------------------------
## 02_run_full.R -- full Monte Carlo, paper setup (parallelised)
##   N = 1000, V = 500, K = 10, P = 8, R = 100 per scenario
##
## Parallelises replicates within each scenario via parallel::mclapply.
## Uses all available physical cores minus 1.
## Progress is logged to stdout and results/sim_full.log.
##
## Outputs:
##   results/full_metrics.rds
##   results/full_metrics.csv
##   results/sim_full.log
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

## ---------------------------------------------------------------------------
## Logging: writes timestamped lines to stdout and to sim_full.log
## ---------------------------------------------------------------------------
LOG_FILE <- file.path(RESULT_DIR, "sim_full.log")
.log_con <- file(LOG_FILE, open = "wt")

.log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(...), "\n")
  cat(msg)
  cat(msg, file = .log_con)
  flush(.log_con)
}

## ---------------------------------------------------------------------------
## Parallel setup
## ---------------------------------------------------------------------------
n_phys  <- max(1L, parallel::detectCores(logical = FALSE))
n_cores <- max(1L, n_phys - 1L)
## mclapply uses fork(), unsupported on Windows -- fall back to serial.
if (.Platform$OS.type == "windows") n_cores <- 1L
.log(sprintf("Cores: %d physical detected, using %d (OS: %s)",
             n_phys, n_cores, .Platform$OS.type))

## ---------------------------------------------------------------------------
## Simulation design
## ---------------------------------------------------------------------------
R_arg <- Sys.getenv("GSCAMM_SIM_R",  unset = "100")
B_arg <- Sys.getenv("GSCAMM_BOOT_B", unset = "200")

DESIGN <- modifyList(DEFAULT_DESIGN, list(
  N = 1000, V = 500, K = 10, P = 8,
  R      = as.integer(R_arg),
  boot_B = as.integer(B_arg)
))
.log("=== Full simulation design ===")
for (ln in capture.output(str(DESIGN))) .log(ln)

## Surface the env-var toggles in the log so each ablation run is
## auto-documented (R/B/C/D/F see README incremental-runs section).
.log(sprintf("=== Toggles ==="))
.log(sprintf("  fit_max_iter   = %d", DESIGN$fit_max_iter))
.log(sprintf("  boot_max_iter  = %d", DESIGN$boot_max_iter))
.log(sprintf("  noise_scale    = %.3f", DESIGN$noise_scale))
.log(sprintf("  use_polish       = %s", DESIGN$use_polish))
.log(sprintf("  sigma2_polish    = %.3f", DESIGN$sigma2_polish))
.log(sprintf("  param_boot_B     = %d (0 = disabled; populates coverage_B_boot_param column)",
             DESIGN$param_boot_B %||% 0L))
.log(sprintf("  use_stm_spectral = %s (counter-factual; STM uses Random init by default)",
             DESIGN$use_stm_spectral))

doc_lengths <- list(baseline       = 1000,
                    high_covariate = 1000,
                    high_sparsity  = 20)

set.seed(2026)
master_seeds <- sample.int(.Machine$integer.max, DESIGN$R)

## ---------------------------------------------------------------------------
## Worker: runs one replicate and returns rows + timing
## (mclapply forks the process so all globals are available)
## ---------------------------------------------------------------------------
.run_worker <- function(r, sc, design, doc_lengths, master_seeds) {
  t0  <- Sys.time()
  sim <- simulate_gscamm(N = design$N, V = design$V,
                         K = design$K, P = design$P,
                         scenario = sc, seed = master_seeds[r],
                         doc_length_mean = doc_lengths[[sc]])
  rows <- tryCatch(
    run_one_replicate(sim, sc, design = design),
    error = function(e) {
      message(sprintf("[%s rep %d] FAILED: %s", sc, r, conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(rows)) {
    rows$replicate <- r
    rows$seed      <- master_seeds[r]
  }
  list(rows = rows, r = r,
       elapsed = as.numeric(Sys.time() - t0, units = "secs"))
}

## ---------------------------------------------------------------------------
## Main loop: scenarios sequential, replicates parallelised in chunks.
## Chunk size = n_cores so we log progress after every wave of workers.
## ---------------------------------------------------------------------------
all_rows   <- list()
t_global   <- Sys.time()
total_done <- 0L
total_reps <- length(DESIGN$scenarios) * DESIGN$R

chunk_ids <- function(R, size)
  split(seq_len(R), ceiling(seq_len(R) / size))

for (sc in DESIGN$scenarios) {
  t_sc <- Sys.time()
  .log(sprintf("========== scenario: %s  [R=%d, %d cores, chunk=%d] ==========",
               sc, DESIGN$R, n_cores, n_cores))

  sc_rows <- list()

  for (chunk in chunk_ids(DESIGN$R, n_cores)) {
    res_chunk <- parallel::mclapply(
      chunk,
      .run_worker,
      sc           = sc,
      design       = DESIGN,
      doc_lengths  = doc_lengths,
      master_seeds = master_seeds,
      mc.cores        = n_cores,
      mc.preschedule  = FALSE   ## better load balancing for variable-cost jobs
    )

    for (res in res_chunk) {
      total_done <- total_done + 1L
      pct        <- 100 * total_done / total_reps
      elapsed_g  <- as.numeric(Sys.time() - t_global, units = "mins")
      .log(sprintf("  %s rep %3d/%d  %5.1fs  |  total: %5.1f min  [%5.1f%%]",
                   sc, res$r, DESIGN$R, res$elapsed, elapsed_g, pct))
      if (!is.null(res$rows))
        sc_rows[[length(sc_rows) + 1L]] <- res$rows
    }

    ## checkpoint after every chunk
    all_rows_now <- c(all_rows, sc_rows)
    if (length(all_rows_now)) {
      saveRDS(do.call(rbind, all_rows_now),
              file.path(RESULT_DIR, "full_metrics_partial.rds"))
    }
  }

  all_rows <- c(all_rows, sc_rows)
  .log(sprintf("  --> scenario %s complete in %.1f min",
               sc, as.numeric(Sys.time() - t_sc, units = "mins")))
}

## ---------------------------------------------------------------------------
## Save final results
## ---------------------------------------------------------------------------
results <- do.call(rbind, all_rows)
## include all method levels actually present so opt-in comparators
## (e.g. stm_spectral) are not coerced to NA by the factor cast.
results$method <- factor(
  results$method,
  levels = intersect(c("gscamm", "gscamm_boot", "lda", "stm", "stm_spectral"),
                     unique(as.character(results$method)))
)

saveRDS(results, file.path(RESULT_DIR, "full_metrics.rds"))
write.csv(results, file.path(RESULT_DIR, "full_metrics.csv"), row.names = FALSE)
unlink(file.path(RESULT_DIR, "full_metrics_partial.rds"))

total_mins <- as.numeric(Sys.time() - t_global, units = "mins")
.log(sprintf("=== Full run complete in %.1f min ===", total_mins))
.log(paste("Results:", file.path(RESULT_DIR, "full_metrics.rds")))
.log(paste("Log:    ", LOG_FILE))
close(.log_con)
