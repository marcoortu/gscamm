## ---------------------------------------------------------------------------
## 04_make_figures.R -- box-plots of the per-replicate metrics, mirroring
## Figures 1-4 of the paper: RMSE_theta, RMSE_B, coverage_B, perplexity.
##
## Reads:  results/full_metrics.rds
## Writes: figures/{rmse_theta,rmse_B,coverage_B,perplexity}.png
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
R$method <- factor(R$method,
                   levels = c("gscamm", "gscamm_boot", "lda", "stm"),
                   labels = c("GSCA-MM\n(plug-in)",
                              "GSCA-MM\n(boot+noise)",
                              "LDA+ALR", "STM"))
R$scenario <- factor(R$scenario,
                     levels = c("baseline", "high_covariate", "high_sparsity"))

palette <- c("#3B82F6", "#1E40AF", "#10B981", "#EF4444")

draw_box <- function(varname, title, ylab,
                     ref_line = NULL, file = NULL, log_y = FALSE) {
  if (is.null(file)) file <- file.path(FIGURE_DIR, paste0(varname, ".png"))
  png(file, width = 1100, height = 700, res = 130)
  op <- par(mar = c(7, 4.2, 3, 1), las = 2)
  on.exit({ par(op); dev.off() })
  ## per-(scenario,method) values
  R$grp <- interaction(R$scenario, R$method, drop = TRUE, sep = " | ")
  vals <- split(R[[varname]], R$grp)
  cols <- rep(palette, each = 3)[match(levels(R$grp), levels(R$grp))]
  log_arg <- if (log_y) "y" else ""
  boxplot(vals, main = title, ylab = ylab, col = cols, log = log_arg,
          outline = FALSE, las = 2, cex.axis = 0.75)
  if (!is.null(ref_line)) abline(h = ref_line, col = "darkred", lty = 2)
  legend("topright", legend = levels(R$method), fill = palette, bty = "n",
         ncol = 2, cex = 0.85)
}

draw_box("rmse_theta", expression(RMSE[theta]),
         expression(RMSE[theta]))
draw_box("rmse_B",     expression(RMSE[B]),
         expression(RMSE[B]))
draw_box("coverage_B", expression("Coverage of 95% CI for "*B),
         "empirical coverage", ref_line = 0.95)
draw_box("perplexity", "Held-in perplexity", "perplexity", log_y = TRUE)
draw_box("time",       "Per-replicate runtime", "seconds (log)", log_y = TRUE)

cat("Figures saved to:\n")
for (f in c("rmse_theta.png", "rmse_B.png",
            "coverage_B.png", "perplexity.png", "time.png"))
  cat("  ", file.path(FIGURE_DIR, f), "\n")
