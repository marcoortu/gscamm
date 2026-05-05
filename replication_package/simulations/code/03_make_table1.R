## ---------------------------------------------------------------------------
## 03_make_table1.R -- aggregate per-replicate metrics into the Table 1
## reproduction (with the additional bootstrap coverage column).
##
## Reads:  results/full_metrics.rds
## Writes: results/table1.csv, results/table1.tex
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

infile <- file.path(RESULT_DIR, "full_metrics.rds")
if (!file.exists(infile))
  stop("Run 02_run_full.R first to produce ", infile)

R <- readRDS(infile)
## detect available methods (stm_spectral is opt-in counter-factual;
## STM means Random init by default in this replication package)
present_methods <- intersect(
  c("gscamm", "gscamm_boot", "lda", "stm", "stm_spectral"),
  unique(as.character(R$method))
)
R$method   <- factor(R$method, levels = present_methods)
R$scenario <- factor(R$scenario,
                     levels = c("baseline", "high_covariate", "high_sparsity"))

## back-compat: columns added in later iterations may be missing in older RDS
if (!"rmse_theta_map" %in% names(R))         R$rmse_theta_map <- R$rmse_theta
if (!"coverage_B_boot_param" %in% names(R))  R$coverage_B_boot_param <- NA_real_
if (!"width_B_boot_param"    %in% names(R))  R$width_B_boot_param    <- NA_real_

## NaN -> NA when the whole column is missing (so the table prints "--"
## via the LaTeX exporter rather than "NaN")
agg_mean <- function(x) { v <- mean(x, na.rm = TRUE); if (is.nan(v)) NA_real_ else v }
agg_sd   <- function(x) { v <- stats::sd(x, na.rm = TRUE); if (is.na(v)) NA_real_ else v }

means <- aggregate(
  cbind(rmse_theta, rmse_theta_R, rmse_theta_map,
        rmse_phi, rmse_B, coverage_B, width_B,
        coverage_B_boot_param, width_B_boot_param,
        perplexity, time)
    ~ scenario + method, data = R, FUN = agg_mean,
  na.action = na.pass)
sds <- aggregate(
  cbind(rmse_theta, rmse_theta_R, rmse_theta_map,
        rmse_phi, rmse_B, coverage_B, width_B,
        coverage_B_boot_param, width_B_boot_param,
        perplexity, time)
    ~ scenario + method, data = R, FUN = agg_sd,
  na.action = na.pass)

tab <- data.frame(
  condition = means$scenario,
  method    = means$method,
  RMSE_theta_mean     = round(means$rmse_theta,     3),
  RMSE_theta_sd       = round(sds$rmse_theta,       3),
  RMSE_theta_R_mean   = round(means$rmse_theta_R,   3),
  RMSE_theta_map_mean = round(means$rmse_theta_map, 3),
  RMSE_theta_map_sd   = round(sds$rmse_theta_map,   3),
  RMSE_phi_mean       = round(means$rmse_phi,       3),
  RMSE_B_mean         = round(means$rmse_B,         3),
  RMSE_B_sd           = round(sds$rmse_B,           3),
  coverage_B_mean     = round(means$coverage_B,     3),
  ## parametric-bootstrap coverage on the gscamm row only (NA elsewhere).
  ## Re-introduced after the model became deterministic so that the
  ## parametric resampling W^(b) ~ Mult(L, theta_hat'Phi_hat) is now
  ## well-specified (the legacy version was removed in earlier reviews
  ## because it under-covered under the old generative model with epsilon).
  coverage_mean_boot  = round(means$coverage_B_boot_param, 3),
  width_B_mean        = round(means$width_B,        3),
  width_mean_boot     = round(means$width_B_boot_param, 3),
  perplexity_mean     = round(means$perplexity,     3),
  time_sec_mean       = round(means$time,           2),
  time_sec_sd         = round(sds$time,             2),
  stringsAsFactors = FALSE
)

## reorder rows: scenario then method
tab <- tab[order(tab$condition, tab$method), ]
print(tab, row.names = FALSE)

write.csv(tab, file.path(RESULT_DIR, "table1.csv"), row.names = FALSE)

## minimal LaTeX export
fmt  <- function(x) if (is.na(x)) "--" else sprintf("%.3f", x)
fmt2 <- function(x) if (is.na(x)) "--" else sprintf("%.2f", x)
lines <- c(
  "\\begin{tabular}{llrrrrrrrrrr}",
  "\\hline",
  paste0("condition & method & RMSE\\_theta\\_mean & RMSE\\_theta\\_sd",
         " & RMSE\\_theta\\_map\\_mean",
         " & RMSE\\_phi\\_mean & RMSE\\_B\\_mean & RMSE\\_B\\_sd",
         " & coverage\\_B\\_mean & coverage\\_mean\\_boot",
         " & perplexity\\_mean",
         " & time\\_sec\\_mean \\\\"),
  "\\hline"
)
for (i in seq_len(nrow(tab))) {
  lines <- c(lines, sprintf(
    "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
    tab$condition[i], tab$method[i],
    fmt(tab$RMSE_theta_mean[i]), fmt(tab$RMSE_theta_sd[i]),
    fmt(tab$RMSE_theta_map_mean[i]),
    fmt(tab$RMSE_phi_mean[i]),
    fmt(tab$RMSE_B_mean[i]), fmt(tab$RMSE_B_sd[i]),
    fmt(tab$coverage_B_mean[i]),
    fmt(tab$coverage_mean_boot[i]),
    fmt(tab$perplexity_mean[i]),
    fmt2(tab$time_sec_mean[i])
  ))
}
lines <- c(lines, "\\hline", "\\end{tabular}")
writeLines(lines, file.path(RESULT_DIR, "table1.tex"))
cat("\nTable saved to:\n  ",
    file.path(RESULT_DIR, "table1.csv"), "\n  ",
    file.path(RESULT_DIR, "table1.tex"), "\n")
