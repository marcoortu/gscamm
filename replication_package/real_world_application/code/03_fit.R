## ---------------------------------------------------------------------------
## 03_fit.R -- fit the GSCA-MM model at the chosen K (default: minimizer of
## held-out perplexity from 02_search_K.R; override via env var).
##
## Reads:  data/processed_dfm.rds, data/processed_X.rds,
##         results/search_K.csv (if present, to pick K)
## Writes: results/fit.rds                (the gscamm object)
##         results/fit_summary.txt        (human-readable summary)
##
## Override K via env var GSCAMM_APP_K. Defaults to the minimizer of
## held-out perplexity from search_K.csv, falling back to K = 8.
## ---------------------------------------------------------------------------

.find_setup <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).*", args, perl = TRUE))
  here <- if (length(m)) dirname(normalizePath(m[1])) else getwd()
  if (!file.exists(file.path(here, "00_setup.R")))
    here <- file.path(getwd(),
                      "replication_package/real_world_application/code")
  file.path(here, "00_setup.R")
}
source(.find_setup())

W <- as.matrix(readRDS(file.path(DATA_DIR, "processed_dfm.rds")))
X <- readRDS(file.path(DATA_DIR, "processed_X.rds"))

## ---------------------------------------------------------------------------
## Choose K
## ---------------------------------------------------------------------------
override <- Sys.getenv("GSCAMM_APP_K", unset = "")
search_K_path <- file.path(RESULT_DIR, "search_K.csv")
if (nzchar(override)) {
  K_use <- as.integer(override)
  cat(sprintf("K = %d (override via GSCAMM_APP_K)\n", K_use))
} else if (file.exists(search_K_path)) {
  res <- read.csv(search_K_path)
  K_use <- res$K[which.min(res$holdout_perplexity)]
  cat(sprintf("K = %d (minimizer of held-out perplexity in %s)\n",
              K_use, search_K_path))
} else {
  K_use <- 8L
  cat(sprintf("K = %d (default; search_K.csv not found)\n", K_use))
}

## ---------------------------------------------------------------------------
## Fit
## ---------------------------------------------------------------------------
t0 <- Sys.time()
ctl <- gscamm_control(max_iter = 200L, tol = 1e-5, polish_max_iter = 50L)
fit <- fit_gscamm(W, X, K = K_use,
                  link = "logistic_normal",
                  gsca_space = "alr",
                  init_phi = "kmeans",
                  polish = "map",
                  sigma2_polish = 0.25,
                  control = ctl,
                  seed = 2026,
                  verbose = FALSE)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("Fit complete: %d iters, converged = %s, %.1fs\n",
            fit$convergence$iterations,
            fit$convergence$converged, elapsed))
cat(sprintf("Final perplexity (in-sample): %.2f\n",
            tail(fit$convergence$perplexity, 1)))

saveRDS(fit, file.path(RESULT_DIR, "fit.rds"))
cat(sprintf("Fit saved -> %s\n", file.path(RESULT_DIR, "fit.rds")))

## ---------------------------------------------------------------------------
## Human-readable summary
## ---------------------------------------------------------------------------
sm <- capture.output({
  cat("=== GSCA-MM fit ===\n")
  cat(sprintf("N = %d, V = %d, K = %d, P = %d\n",
              fit$N, fit$V, fit$K, fit$P))
  cat(sprintf("Link: %s    GSCA space: %s    init_phi: %s    polish: %s\n",
              fit$link, fit$gsca_space, fit$init_phi, fit$polish))
  cat(sprintf("EM iterations: %d (converged = %s)\n",
              fit$convergence$iterations, fit$convergence$converged))
  cat(sprintf("Final in-sample perplexity: %.2f\n",
              tail(fit$convergence$perplexity, 1)))
  cat(sprintf("MAP polish: %d/%d rows converged\n",
              fit$convergence$polish$n_converged %||% NA, fit$N))
  cat("\nsigma2_k (estimator from R-residuals, diagnostic):\n")
  print(round(fit$sigma2_k, 4))
  cat("\nsigma2_polish (used in MAP):\n")
  print(round(fit$sigma2_polish, 4))
  cat("\nB_minus (P x (K-1)):\n")
  print(round(fit$B_minus, 3))
})
writeLines(sm, file.path(RESULT_DIR, "fit_summary.txt"))
cat(sprintf("Summary saved -> %s\n",
            file.path(RESULT_DIR, "fit_summary.txt")))
