## ---------------------------------------------------------------------------
## 05_run_ablation.R -- drives the 4 incremental ablation runs (B, C, D, F)
## end-to-end with a single Rscript invocation.
##
## Each child process re-sources 00_setup.R and runs the full simulation
## under its own env-var configuration. Between runs, this driver renames
## the result files with a run-specific suffix so all outputs survive on
## disk for differential analysis. If a `full_metrics.rds` already exists
## at the start, it is snapshotted to `full_metrics_runA.rds` first so the
## pre-ablation baseline is preserved.
##
## Outputs (per run X in {A, B, C, D, F}):
##   results/full_metrics_runX.rds
##   results/full_metrics_runX.csv
##   results/sim_full_runX.log
## Plus:
##   results/ablation_summary.csv  (per-run status + elapsed minutes)
##
## Usage (project root, server):
##   R CMD INSTALL .                                   # install package first!
##   Rscript replication_package/simulations/code/05_run_ablation.R
##
## NOTE on installation: each ablation run is a fresh Rscript subprocess
## that loads gscamm via library(), so any source-only changes (devtools::
## load_all) are NOT visible to the children. The package MUST be installed
## (R CMD INSTALL . or devtools::install()) before this driver is launched.
##
## Optional env-var overrides (also forwarded to child runs):
##   GSCAMM_SIM_R   replicates per scenario (default 100)
##   GSCAMM_BOOT_B  bootstrap reps         (default 200)
##
## Approximate cost on linux02-class machine at R=100, B=200:
##   runB ~30 min, runC ~50 min, runD ~55 min, runF ~50 min  (total ~3 h)
## ---------------------------------------------------------------------------

.find_self <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).*", args, perl = TRUE))
  if (length(m)) {
    normalizePath(m[1])
  } else if (!is.null(sys.frames()) && length(sys.frames())) {
    of <- sys.frames()[[1]]$ofile
    if (!is.null(of)) normalizePath(of) else NA_character_
  } else {
    NA_character_
  }
}

HERE <- tryCatch(dirname(.find_self()), error = function(e) NA_character_)
if (is.na(HERE) || !dir.exists(HERE))
  HERE <- file.path(getwd(), "replication_package/simulations/code")

REPL_ROOT  <- normalizePath(file.path(HERE, ".."))
RESULT_DIR <- file.path(REPL_ROOT, "results")
RUN_FULL   <- file.path(HERE, "02_run_full.R")
if (!file.exists(RUN_FULL))
  stop("Cannot locate 02_run_full.R at: ", RUN_FULL)
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

## ---------------------------------------------------------------------------
## Configuration of the 4 incremental ablation runs (B/C/D/F).
## See replication_package/simulations/README.md for the rationale of each.
## ---------------------------------------------------------------------------
## NOTE: this driver was used to validate the design changes that became
## the new defaults (MAP polish on, STM uses Random init, noise_scale=1.0,
## fit_max_iter=200). After validation, STM Random became the canonical
## comparator and Spectral is the opt-in counterfactual; the env-var was
## renamed accordingly. The original ablation logs and result files are
## preserved on disk (full_metrics_run{B,C,D,F}.{rds,csv}). The configs
## below replay the same ablation under the new naming for completeness;
## for a production run prefer Rscript 02_run_full.R directly.
ABLATION_CONFIGS <- list(
  list(tag  = "runB",
       desc = "fit_max_iter=200, no polish, no stm_spectral",
       env  = list(GSCAMM_USE_POLISH       = "0",
                   GSCAMM_USE_STM_SPECTRAL = "0",
                   GSCAMM_NOISE_SCALE      = "2.0",
                   GSCAMM_FIT_MAX_ITER     = "200",
                   GSCAMM_BOOT_MAX_ITER    = "60")),
  list(tag  = "runC",
       desc = "+stm_spectral comparator (isolate Spectral-init contribution)",
       env  = list(GSCAMM_USE_POLISH       = "0",
                   GSCAMM_USE_STM_SPECTRAL = "1",
                   GSCAMM_NOISE_SCALE      = "2.0",
                   GSCAMM_FIT_MAX_ITER     = "200",
                   GSCAMM_BOOT_MAX_ITER    = "60")),
  list(tag  = "runD",
       desc = "+MAP polish (apples-to-apples theta estimator)",
       env  = list(GSCAMM_USE_POLISH       = "1",
                   GSCAMM_USE_STM_SPECTRAL = "1",
                   GSCAMM_NOISE_SCALE      = "2.0",
                   GSCAMM_FIT_MAX_ITER     = "200",
                   GSCAMM_BOOT_MAX_ITER    = "60",
                   GSCAMM_SIGMA2_POLISH    = "0.25")),
  list(tag  = "runF",
       desc = "+noise_scale=1.0 (recalibrate boot coverage)",
       env  = list(GSCAMM_USE_POLISH       = "1",
                   GSCAMM_USE_STM_SPECTRAL = "1",
                   GSCAMM_NOISE_SCALE      = "1.0",
                   GSCAMM_FIT_MAX_ITER     = "200",
                   GSCAMM_BOOT_MAX_ITER    = "60",
                   GSCAMM_SIGMA2_POLISH    = "0.25"))
)

R_arg <- Sys.getenv("GSCAMM_SIM_R",  unset = "100")
B_arg <- Sys.getenv("GSCAMM_BOOT_B", unset = "200")

## ---------------------------------------------------------------------------
## Snapshot the existing results/full_metrics.{rds,csv} (if any) as the
## pre-ablation baseline (runA). Skip if a runA snapshot already exists --
## never overwrite an earlier baseline.
## ---------------------------------------------------------------------------
.snapshot_if_present <- function(src_name, dst_name) {
  src <- file.path(RESULT_DIR, src_name)
  dst <- file.path(RESULT_DIR, dst_name)
  if (file.exists(src) && !file.exists(dst)) {
    if (file.copy(src, dst))
      cat(sprintf("Baseline snapshot: %s -> %s\n", src_name, dst_name))
  }
}
.snapshot_if_present("full_metrics.rds", "full_metrics_runA.rds")
.snapshot_if_present("full_metrics.csv", "full_metrics_runA.csv")
.snapshot_if_present("sim_full.log",     "sim_full_runA.log")

## ---------------------------------------------------------------------------
## Driver loop
## ---------------------------------------------------------------------------
.banner <- function(s) cat(sprintf("\n%s\n%s\n%s\n", strrep("=", 72), s,
                                   strrep("=", 72)))

started <- Sys.time()
.banner(sprintf("[%s] Ablation start: %d runs  R=%s  B=%s",
                format(started, "%Y-%m-%d %H:%M:%S"),
                length(ABLATION_CONFIGS), R_arg, B_arg))

results <- data.frame(tag = character(), status = character(),
                      elapsed_min = numeric(),
                      stringsAsFactors = FALSE)

## stale-file names that 02_run_full.R writes to (cleared before each run)
STALE_FILES <- c("full_metrics.rds", "full_metrics.csv",
                 "full_metrics_partial.rds", "sim_full.log")

for (cfg in ABLATION_CONFIGS) {
  .banner(sprintf("%s : %s", cfg$tag, cfg$desc))

  ## clear any previous-run leftover so we don't mistake stale state for new
  for (fn in STALE_FILES) {
    p <- file.path(RESULT_DIR, fn)
    if (file.exists(p)) file.remove(p)
  }
  ## also clear any previously-tagged outputs FOR THIS RUN so a prior
  ## attempt's leftover cannot be mistaken for the new attempt's result
  ## (runA is preserved -- the pre-ablation baseline -- and is never touched)
  for (fn in c(sprintf("full_metrics_%s.rds", cfg$tag),
               sprintf("full_metrics_%s.csv", cfg$tag),
               sprintf("sim_full_%s.log",     cfg$tag),
               sprintf("full_metrics_%s_partial.rds", cfg$tag))) {
    p <- file.path(RESULT_DIR, fn)
    if (file.exists(p)) file.remove(p)
  }

  ## env vars: set in the parent process before spawning the child (child
  ## inherits them). This is portable; system2(env=) is unreliable on
  ## Windows (treated as positional args). All keys touched here are
  ## unset at end of iteration so a later run cannot inherit a stale toggle.
  child_env <- c(unlist(cfg$env),
                 GSCAMM_SIM_R  = R_arg,
                 GSCAMM_BOOT_B = B_arg)
  cat("env: ", paste(sprintf("%s=%s", names(child_env), child_env),
                     collapse = "  "), "\n", sep = "")
  do.call(Sys.setenv, as.list(child_env))

  rscript_bin <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows") rscript_bin <- paste0(rscript_bin, ".exe")

  t0 <- Sys.time()
  rc <- tryCatch(
    system2(rscript_bin,
            args = shQuote(RUN_FULL),
            stdout = "", stderr = ""),
    error = function(e) {
      message("** ", cfg$tag, " threw: ", conditionMessage(e)); -1L
    }
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "mins")

  ## unset the toggles we set so a subsequent run starts from defaults
  Sys.unsetenv(names(child_env))

  status <- if (identical(as.integer(rc), 0L)) "ok"
            else sprintf("failed(rc=%s)", as.character(rc))
  cat(sprintf("\n[%s] %s -> %s in %.1f min\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              cfg$tag, status, elapsed))

  ## rename whatever survives so each run yields a tagged artifact even on
  ## partial failure (the partial RDS is useful for a post-mortem)
  rename_pairs <- list(
    c("full_metrics.rds",         sprintf("full_metrics_%s.rds", cfg$tag)),
    c("full_metrics.csv",         sprintf("full_metrics_%s.csv", cfg$tag)),
    c("sim_full.log",             sprintf("sim_full_%s.log",     cfg$tag)),
    c("full_metrics_partial.rds", sprintf("full_metrics_%s_partial.rds",
                                          cfg$tag))
  )
  for (pair in rename_pairs) {
    src <- file.path(RESULT_DIR, pair[1])
    dst <- file.path(RESULT_DIR, pair[2])
    if (file.exists(src)) {
      ## file.rename can fail across volumes; fall back to copy + remove
      ok <- tryCatch(file.rename(src, dst), warning = function(w) FALSE)
      if (!isTRUE(ok)) {
        ok <- file.copy(src, dst, overwrite = TRUE) && file.remove(src)
      }
      if (isTRUE(ok)) cat(sprintf("  saved -> %s\n", basename(dst)))
    }
  }

  results <- rbind(results,
                   data.frame(tag = cfg$tag, status = status,
                              elapsed_min = round(elapsed, 1),
                              stringsAsFactors = FALSE))

  ## persist the summary after every run so a crash leaves a partial trail
  write.csv(results, file.path(RESULT_DIR, "ablation_summary.csv"),
            row.names = FALSE)
}

total_min <- as.numeric(Sys.time() - started, units = "mins")
.banner(sprintf("Ablation complete in %.1f min", total_min))
print(results, row.names = FALSE)

cat(sprintf("\nArtifacts in: %s\n", RESULT_DIR))
cat("  full_metrics_run{A,B,C,D,F}.{rds,csv}\n")
cat("  sim_full_run{A,B,C,D,F}.log\n")
cat("  ablation_summary.csv\n")
