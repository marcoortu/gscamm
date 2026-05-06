## ---------------------------------------------------------------------------
## 05_covariate_effects.R -- estimate and visualise the structural covariate
## effects on topic prevalence.
##
## Produces:
##   figures/effects_forest.png      forest plot, faceted by covariate
##   figures/effects_heatmap.png     covariates x topics, with topic
##                                   hierarchical clustering on the columns
##   figures/prevalence_by_year.png  predicted topic prevalence by year
##                                   (party-marginal & candidate-marginal)
##   results/effects_long.csv        long-format coefficient table
##
## Reads:  results/fit.rds, data/processed_X.rds, data/processed_df.rds
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

suppressPackageStartupMessages({
  library(scales)
  library(ggrepel)
})

## Inlined tidytext helpers: reorder discrete factors within facets
## without taking a tidytext dependency.
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

fit <- readRDS(file.path(RESULT_DIR, "fit.rds"))
X   <- readRDS(file.path(DATA_DIR, "processed_X.rds"))
df  <- readRDS(file.path(DATA_DIR, "processed_df.rds"))

## ---------------------------------------------------------------------------
## Compute plug-in covariate effects (uses MAP polish in fit$Theta_map)
## ---------------------------------------------------------------------------
eff <- covariate_effects(fit, ref = fit$K, level = 0.95, adjust = "BH",
                         standardize = TRUE,
                         use = "theta", var_method = "hybrid")
co <- eff$coefficients
co <- co[co$covariate != "(Intercept)", , drop = FALSE]

## human-friendly covariate labels for the figures
.cov_label <- c(
  "year2020"                       = "Year = 2020 (vs 2019)",
  "year2021"                       = "Year = 2021 (vs 2019)",
  "partydemocratic"                = "Party = Democratic (vs Republican)",
  "candidate_typevice_president"   = "Candidate = VP (vs President)"
)
co$covariate_label <- factor(.cov_label[co$covariate],
                             levels = unname(.cov_label))
co$component_label <- factor(co$component,
                             levels = paste0("comp", seq_len(fit$K - 1L)),
                             labels = paste0("Topic ",
                                             seq_len(fit$K - 1L)))
co$significant <- co$p.adj < 0.05
write.csv(co, file.path(RESULT_DIR, "effects_long.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Figure 1: forest plot, faceted by covariate.
## Within each facet, topics are sorted by effect size for visual ranking.
## Significant effects (BH-adjusted p < 0.05) are colored, others greyed out.
## ---------------------------------------------------------------------------
co_plot <- co %>%
  group_by(covariate_label) %>%
  mutate(.order = rank(estimate, ties.method = "first")) %>%
  ungroup() %>%
  mutate(component_ord = factor(
    paste(covariate_label, component_label, sep = "##"),
    levels = unique(paste(covariate_label, component_label,
                          sep = "##")[order(covariate_label, .order)])
  ))

p_forest <- ggplot(co_plot,
                   aes(x = estimate,
                       y = reorder_within(component_label,
                                          estimate, covariate_label),
                       colour = ifelse(significant,
                                       ifelse(estimate > 0, "pos", "neg"),
                                       "ns"))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25,
                 linewidth = 0.5) +
  geom_point(size = 2.2) +
  facet_wrap(~ covariate_label, scales = "free", ncol = 2) +
  scale_y_reordered() +
  scale_colour_manual(values = c("pos" = "#1E40AF",
                                  "neg" = "#B91C1C",
                                  "ns"  = "grey60"),
                      labels = c("pos" = "Positive (BH < 0.05)",
                                 "neg" = "Negative (BH < 0.05)",
                                 "ns"  = "Non-significant"),
                      breaks = c("pos", "neg", "ns"),
                      name = NULL) +
  labs(title = "Covariate effects on topic prevalence (ALR scale)",
       subtitle = "Plug-in WLS with MAP-polish theta. Bars are 95% CIs.",
       x = "ALR coefficient (log-odds vs reference topic)",
       y = NULL) +
  theme_paper()

## reorder_within / scale_y_reordered helpers (avoid pulling in tidytext just for this)
ggsave(file.path(FIGURE_DIR, "effects_forest.png"),
       plot = p_forest,
       width = 10, height = 1.5 + 0.45 * (fit$K - 1) * 2,
       dpi = 150, limitsize = FALSE)
cat(sprintf("Forest plot -> %s\n",
            file.path(FIGURE_DIR, "effects_forest.png")))

## ---------------------------------------------------------------------------
## Figure 2: heatmap of covariate x topic effects, with hierarchical
## clustering on topics so that topics with similar effect signatures
## are adjacent. Cell colour = ALR coefficient (log-odds), with 0 mapped to
## white via a diverging palette.
## ---------------------------------------------------------------------------
mat <- matrix(NA_real_,
              nrow = length(.cov_label), ncol = fit$K - 1L,
              dimnames = list(unname(.cov_label),
                              paste0("Topic ", seq_len(fit$K - 1L))))
for (i in seq_len(nrow(co))) {
  r <- match(co$covariate_label[i], rownames(mat))
  c <- match(co$component_label[i], colnames(mat))
  if (!is.na(r) && !is.na(c)) mat[r, c] <- co$estimate[i]
}

## hierarchical cluster on topics (columns)
ord_cols <- if (ncol(mat) >= 2) {
  d <- dist(t(mat))
  hc <- hclust(d, method = "ward.D2")
  hc$order
} else seq_len(ncol(mat))

mat_ord <- mat[, ord_cols, drop = FALSE]
heat_df <- expand.grid(covariate = rownames(mat_ord),
                       topic = colnames(mat_ord),
                       stringsAsFactors = FALSE) %>%
  mutate(value = as.vector(mat_ord),
         topic = factor(topic, levels = colnames(mat_ord)),
         covariate = factor(covariate, levels = rev(rownames(mat_ord))))

vmax <- max(abs(heat_df$value), na.rm = TRUE)
p_heat <- ggplot(heat_df, aes(x = topic, y = covariate, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", value)),
            size = 3, colour = "black") +
  scale_fill_gradient2(low = "#B91C1C", mid = "white", high = "#1E40AF",
                       midpoint = 0,
                       limits = c(-vmax, vmax),
                       name = "ALR coef") +
  labs(title = "Covariate x Topic effect heatmap",
       subtitle = "Topics ordered by hierarchical clustering of effect signatures",
       x = NULL, y = NULL) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_blank())

ggsave(file.path(FIGURE_DIR, "effects_heatmap.png"),
       plot = p_heat,
       width = max(8, 0.7 * fit$K + 4), height = 4.5, dpi = 150)
cat(sprintf("Heatmap -> %s\n",
            file.path(FIGURE_DIR, "effects_heatmap.png")))

## ---------------------------------------------------------------------------
## Figure 3: predicted topic prevalence over years, marginalised over the
## empirical distribution of party x candidate_type. For each year level
## we set X[, year_*] to the appropriate dummy values and average theta_i
## across the documents (which keeps the empirical mix of party / VP).
## ---------------------------------------------------------------------------
years <- c("2019", "2020", "2021")
X_means_by_year <- vapply(years, function(yr) {
  Xn <- X
  Xn[, "year2020"] <- as.integer(yr == "2020")
  Xn[, "year2021"] <- as.integer(yr == "2021")
  ## standardize the same way the fit did
  ctr <- attr(scale(X), "scaled:center")
  scl <- attr(scale(X), "scaled:scale")
  scl[scl == 0] <- 1
  Xn_std <- sweep(sweep(Xn, 2, ctr, "-"), 2, scl, "/")
  ## structural eta = X_std %*% B_minus
  eta <- Xn_std %*% fit$B_minus
  eta_full <- cbind(eta, 0)
  m <- apply(eta_full, 1, max)
  th <- exp(eta_full - m); th <- th / rowSums(th)
  colMeans(th)
}, numeric(fit$K))

prev_df <- as.data.frame(t(X_means_by_year))
colnames(prev_df) <- paste0("Topic ", seq_len(fit$K))
prev_df$year <- factor(rownames(prev_df), levels = years)
prev_long <- prev_df %>%
  tidyr::pivot_longer(cols = starts_with("Topic"),
                      names_to = "topic", values_to = "prevalence")

p_prev <- ggplot(prev_long, aes(x = year, y = prevalence,
                                colour = topic, group = topic)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_colour_viridis_d(option = "D", begin = 0.1, end = 0.9) +
  labs(title = "Predicted mean topic prevalence by year",
       subtitle = "Marginalised over the empirical party x candidate-type mix",
       x = NULL, y = "Mean theta",
       colour = NULL) +
  theme_paper()

ggsave(file.path(FIGURE_DIR, "prevalence_by_year.png"),
       plot = p_prev, width = 7.5, height = 4.5, dpi = 150)
cat(sprintf("Prevalence-by-year -> %s\n",
            file.path(FIGURE_DIR, "prevalence_by_year.png")))

cat(sprintf("Effects table -> %s\n",
            file.path(RESULT_DIR, "effects_long.csv")))
