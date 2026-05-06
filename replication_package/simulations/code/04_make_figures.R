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

## back-compat: rmse_theta_map missing in older RDS files
if (!"rmse_theta_map" %in% names(R)) R$rmse_theta_map <- R$rmse_theta

## Drop gscamm_boot rows globally: the figures present GSCA-MM as a single
## method (the plug-in fit). The bootstrap row's RMSE/perplexity values are
## identical to the plug-in (it inherits Theta/Phi), and its coverage/width
## differences live in the table only.
R <- R[as.character(R$method) != "gscamm_boot", , drop = FALSE]

## detect available methods (stm_spectral is opt-in counter-factual;
## STM means Random init by default in this replication package)
present_methods <- intersect(
  c("gscamm", "lda", "stm", "stm_spectral"),
  unique(as.character(R$method))
)
method_labels <- c(gscamm       = "GSCA-MM",
                   lda          = "LDA+ALR",
                   stm          = "STM",
                   stm_spectral = "STM\n(Spectral, warm)")
R$method   <- factor(R$method, levels = present_methods,
                     labels = method_labels[present_methods])
R$scenario <- factor(R$scenario,
                     levels = c("baseline", "high_covariate", "high_sparsity"))

## High-contrast palette: GSCA-MM in deep navy; STM (Random) in deep red
## with the optional Spectral warm-start variant in light red so the
## warm-start contribution is immediately visible.
PALETTE <- c("GSCA-MM"                = "#1E3A8A",
             "LDA+ALR"                = "#10B981",
             "STM"                    = "#EF4444",
             "STM\n(Spectral, warm)"  = "#FCA5A5")

draw_box <- function(varname, title, ylab,
                     ref_line = NULL, file = NULL, log_y = FALSE) {
  if (is.null(file)) file <- file.path(FIGURE_DIR, paste0(varname, ".png"))
  d <- R
  d$method   <- droplevels(d$method)
  d$scenario <- droplevels(d$scenario)

  png(file, width = 1200, height = 900, res = 130)
  op <- par(mar = c(9.0, 4.6, 5.0, 1.5), xpd = FALSE)
  on.exit({ par(op); dev.off() })

  d$grp <- interaction(d$scenario, d$method, drop = TRUE, sep = " | ")
  vals  <- split(d[[varname]], d$grp)
  pal_used <- PALETTE[levels(d$method)]
  cols <- rep(pal_used, each = nlevels(d$scenario))

  log_arg <- if (log_y) "y" else ""
  ## suppress default x-axis labels; we draw them at 45 degrees below
  boxplot(vals, main = "", ylab = ylab, col = cols, log = log_arg,
          outline = FALSE, xaxt = "n", las = 1, cex.axis = 0.85)
  if (!is.null(ref_line)) abline(h = ref_line, col = "darkred", lty = 2)

  ## 45-degree rotated x-axis labels (flatten internal newlines)
  x_labels <- gsub("\n", " ", names(vals))
  axis(1, at = seq_along(vals), labels = FALSE, tick = TRUE)
  usr <- par("usr")
  y_pos <- if (log_y) 10^(usr[3] - 0.04 * (usr[4] - usr[3])) else
                       usr[3] - 0.03 * (usr[4] - usr[3])
  text(x = seq_along(vals), y = y_pos, labels = x_labels,
       srt = 45, adj = c(1, 1), xpd = TRUE, cex = 0.78)

  ## title and legend both placed in the (enlarged) top margin so they do
  ## not occlude any boxplot
  title(main = title, line = 3.5, cex.main = 1.2)
  par(xpd = NA)
  legend_labels <- gsub("\n", " ", levels(d$method))
  legend("top", inset = c(0, -0.10),
         legend = legend_labels, fill = pal_used,
         horiz = TRUE, bty = "n", cex = 0.85)
}

draw_box("rmse_theta", expression(RMSE[theta]*" (structural / posterior)"),
         expression(RMSE[theta]))
draw_box("rmse_theta_map",
         expression(RMSE[theta]),
         expression(RMSE[theta]))
draw_box("rmse_B",     expression(RMSE[B]),
         expression(RMSE[B]))
draw_box("coverage_B", expression("Coverage of 95% CI for "*B),
         "empirical coverage", ref_line = 0.95)
draw_box("perplexity", "Held-in perplexity", "perplexity",
         log_y = TRUE)
draw_box("time",       "Per-replicate runtime", "seconds (log)", log_y = TRUE)

cat("Figures saved to:\n")
for (f in c("rmse_theta.png", "rmse_theta_map.png", "rmse_B.png",
            "coverage_B.png", "perplexity.png", "time.png"))
  cat("  ", file.path(FIGURE_DIR, f), "\n")
