## ---------------------------------------------------------------------------
## 02_search_K.R -- pick the number of latent components by held-out
## perplexity. Fits the GSCA-MM model for each K in a grid and reports
## both in-sample and held-out perplexity (using a 10% token holdout).
##
## Reads:  data/processed_dfm.rds, data/processed_X.rds
## Writes: results/search_K.csv      (table of K vs perplexity)
##         figures/search_K.png      (plot)
##
## Override the search range and holdout via env vars:
##   GSCAMM_KMIN, GSCAMM_KMAX (default 2 ... 18)
##   GSCAMM_HOLDOUT (default 0.10)
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
cat(sprintf("Loaded DFM %d x %d, X %d x %d\n",
            nrow(W), ncol(W), nrow(X), ncol(X)))

K_min   <- as.integer(Sys.getenv("GSCAMM_KMIN",   unset = "2"))
K_max   <- as.integer(Sys.getenv("GSCAMM_KMAX",   unset = "18"))
holdout <- as.numeric(Sys.getenv("GSCAMM_HOLDOUT", unset = "0.10"))
Ks <- seq(K_min, K_max)

t0 <- Sys.time()
cat(sprintf("Searching K in [%d, %d], holdout = %.2f, fit_max_iter = 100\n",
            K_min, K_max, holdout))
ctl <- gscamm_control(max_iter = 100L, tol = 1e-4)
res <- search_optimal_components(W, X, Ks = Ks,
                                 link = "logistic_normal",
                                 control = ctl,
                                 holdout = holdout,
                                 seed = 2026,
                                 verbose = TRUE)
elapsed <- as.numeric(Sys.time() - t0, units = "mins")
cat(sprintf("Done in %.1f min\n", elapsed))
print(res)

write.csv(res, file.path(RESULT_DIR, "search_K.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Plot: in-sample vs held-out perplexity, with a marker on the K we'll use
## ---------------------------------------------------------------------------
res_long <- tidyr::pivot_longer(res,
  cols = c("perplexity", "holdout_perplexity"),
  names_to = "kind", values_to = "value") %>%
  mutate(kind = factor(kind,
                       levels = c("perplexity", "holdout_perplexity"),
                       labels = c("In-sample", "Held-out (10%)")))

## Pick K* = K minimizing held-out perplexity (or the elbow if ties)
K_star <- res$K[which.min(res$holdout_perplexity)]
cat(sprintf("\nMinimum held-out perplexity at K = %d\n", K_star))

p <- ggplot(res_long, aes(x = K, y = value, colour = kind, group = kind)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_vline(xintercept = K_star, linetype = "dashed", colour = "grey40") +
  annotate("text", x = K_star, y = max(res_long$value, na.rm = TRUE),
           label = sprintf("K* = %d", K_star), hjust = -0.1, vjust = 1,
           colour = "grey20") +
  scale_colour_manual(values = c("In-sample" = "#94A3B8",
                                  "Held-out (10%)" = "#1E3A8A")) +
  scale_x_continuous(breaks = Ks) +
  labs(title = "Choosing the number of components K",
       subtitle = "Held-out perplexity selects the operating point",
       x = "Number of components K",
       y = "Perplexity",
       colour = NULL) +
  theme_paper()

ggsave(file.path(FIGURE_DIR, "search_K.png"),
       plot = p, width = 7.5, height = 4.5, dpi = 150)
cat(sprintf("Plot saved -> %s\n", file.path(FIGURE_DIR, "search_K.png")))
cat(sprintf("Table saved -> %s\n", file.path(RESULT_DIR, "search_K.csv")))
