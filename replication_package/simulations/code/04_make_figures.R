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

## High-contrast palette: the two GSCA-MM variants use distinct hues
## (light blue / deep navy) rather than two adjacent blues, so the
## plug-in and bootstrap boxes are visually separable when both are shown.
PALETTE <- c("GSCA-MM\n(plug-in)"    = "#60A5FA",
             "GSCA-MM\n(boot+noise)" = "#1E3A8A",
             "LDA+ALR"               = "#10B981",
             "STM"                   = "#EF4444")

draw_box <- function(varname, title, ylab,
                     ref_line = NULL, file = NULL, log_y = FALSE,
                     drop_boot = FALSE) {
  if (is.null(file)) file <- file.path(FIGURE_DIR, paste0(varname, ".png"))
  d <- R
  ## For metrics where the bootstrap wrapper inherits the plug-in fit's
  ## Theta and Phi (rmse_theta, rmse_phi, perplexity), the boot row is a
  ## duplicate of the plug-in row. Drop it so the boxes don't overlap.
  if (drop_boot) {
    d <- d[d$method != "GSCA-MM\n(boot+noise)", , drop = FALSE]
  }
  d$method   <- droplevels(d$method)
  d$scenario <- droplevels(d$scenario)

  png(file, width = 1200, height = 780, res = 130)
  op <- par(mar = c(7.5, 4.6, 5.0, 1.5), las = 2, xpd = FALSE)
  on.exit({ par(op); dev.off() })

  d$grp <- interaction(d$scenario, d$method, drop = TRUE, sep = " | ")
  vals  <- split(d[[varname]], d$grp)
  pal_used <- PALETTE[levels(d$method)]
  cols <- rep(pal_used, each = nlevels(d$scenario))

  log_arg <- if (log_y) "y" else ""
  ## title is drawn separately so we can leave room for an out-of-plot legend
  boxplot(vals, main = "", ylab = ylab, col = cols, log = log_arg,
          outline = FALSE, las = 2, cex.axis = 0.78)
  if (!is.null(ref_line)) abline(h = ref_line, col = "darkred", lty = 2)

  ## title and legend both placed in the (enlarged) top margin so they do
  ## not occlude any boxplot
  title(main = title, line = 3.5, cex.main = 1.2)
  par(xpd = NA)
  legend_labels <- gsub("\n", " ", levels(d$method))
  legend("top", inset = c(0, -0.10),
         legend = legend_labels, fill = pal_used,
         horiz = TRUE, bty = "n", cex = 0.85)
}

draw_box("rmse_theta", expression(RMSE[theta]),
         expression(RMSE[theta]), drop_boot = TRUE)
draw_box("rmse_B",     expression(RMSE[B]),
         expression(RMSE[B]))
draw_box("coverage_B", expression("Coverage of 95% CI for "*B),
         "empirical coverage", ref_line = 0.95)
draw_box("perplexity", "Held-in perplexity", "perplexity",
         log_y = TRUE, drop_boot = TRUE)
draw_box("time",       "Per-replicate runtime", "seconds (log)", log_y = TRUE)

cat("Figures saved to:\n")
for (f in c("rmse_theta.png", "rmse_B.png",
            "coverage_B.png", "perplexity.png", "time.png"))
  cat("  ", file.path(FIGURE_DIR, f), "\n")
