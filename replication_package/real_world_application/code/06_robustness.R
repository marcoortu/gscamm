## ---------------------------------------------------------------------------
## 06_robustness.R -- robustness check via LDA+ALR and STM-Random.
##
## Re-fits LDA and STM on the same DFM/covariates and applies the same
## ALR-WLS covariate-effects estimator. Compares the three methods through:
##   (a) per-(covariate, topic) coefficient comparison plot, paired by
##       topic alignment (Hungarian on Phi);
##   (b) per-method count of significant effects per covariate.
##
## Model-to-model fairness: STM uses init.type = "Random" by default so the
## anchor-words spectral warm-start does not give STM an initialisation
## advantage. (See note in simulations/README.md.)
##
## Reads:  results/fit.rds, data/processed_dfm.rds, data/processed_X.rds,
##         data/processed_df.rds
## Writes: figures/method_comparison_effects.png
##         figures/method_comparison_signif_count.png
##         results/effects_lda.csv, results/effects_stm.csv,
##         results/effects_all.csv
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
  library(topicmodels)
  library(stm)
  library(slam)
})

fit <- readRDS(file.path(RESULT_DIR, "fit.rds"))
W   <- as.matrix(readRDS(file.path(DATA_DIR, "processed_dfm.rds")))
X   <- readRDS(file.path(DATA_DIR, "processed_X.rds"))
df  <- readRDS(file.path(DATA_DIR, "processed_df.rds"))
K   <- fit$K

## standardize X the same way fit_gscamm did
ctr <- colMeans(X)
scl <- apply(X, 2, sd); scl[scl == 0] <- 1
X_std <- sweep(sweep(X, 2, ctr, "-"), 2, scl, "/")

## ---------------------------------------------------------------------------
## 1. Fit LDA on the integer DFM (round to int defensively)
## ---------------------------------------------------------------------------
cat("Fitting LDA (k = ", K, ")...\n", sep = "")
W_int <- round(W); storage.mode(W_int) <- "integer"
dtm_lda <- slam::as.simple_triplet_matrix(W_int)
t0 <- Sys.time()
lda_fit <- topicmodels::LDA(dtm_lda, k = K, method = "VEM",
                            control = list(alpha = 0.1, seed = 2026))
cat(sprintf("  LDA fit in %.1fs\n",
            as.numeric(Sys.time() - t0, units = "secs")))
Phi_lda   <- exp(lda_fit@beta); Phi_lda <- Phi_lda / rowSums(Phi_lda)
Theta_lda <- lda_fit@gamma

## ---------------------------------------------------------------------------
## 2. Fit STM with Random init (canonical comparator) and prevalence on
##    the same covariates
## ---------------------------------------------------------------------------
cat("Fitting STM (init.type = Random, k = ", K, ")...\n", sep = "")
.dtm_to_stm <- function(W_int) {
  V <- ncol(W_int); vocab <- colnames(W_int)
  if (is.null(vocab)) vocab <- sprintf("v%d", seq_len(V))
  documents <- lapply(seq_len(nrow(W_int)), function(i) {
    nz <- which(W_int[i, ] > 0)
    if (!length(nz)) return(matrix(0L, nrow = 2L, ncol = 0L))
    rbind(as.integer(nz), as.integer(W_int[i, nz]))
  })
  list(documents = documents, vocab = vocab)
}
stm_in <- .dtm_to_stm(W_int)
proc <- stm::prepDocuments(stm_in$documents, stm_in$vocab,
                           lower.thresh = 0L, verbose = FALSE)
kept_idx <- match(proc$vocab, colnames(W_int))

X_df <- as.data.frame(X)
prev_form <- as.formula(paste("~", paste(colnames(X_df), collapse = "+")))
t0 <- Sys.time()
set.seed(2026)
stm_fit <- stm::stm(documents = proc$documents, vocab = proc$vocab,
                    K = K, prevalence = prev_form, data = X_df,
                    max.em.its = 75, init.type = "Random", verbose = FALSE)
cat(sprintf("  STM fit in %.1fs\n",
            as.numeric(Sys.time() - t0, units = "secs")))
Theta_stm <- stm_fit$theta
Phi_stm_kept <- exp(stm_fit$beta$logbeta[[1]])
Phi_stm_kept <- Phi_stm_kept / rowSums(Phi_stm_kept)
Phi_stm <- matrix(0, nrow = K, ncol = ncol(W))
Phi_stm[, kept_idx] <- Phi_stm_kept
rs <- rowSums(Phi_stm); rs[rs == 0] <- 1; Phi_stm <- Phi_stm / rs

## ---------------------------------------------------------------------------
## 3. Hungarian alignment of LDA / STM topics to GSCA-MM topics
## ---------------------------------------------------------------------------
perm_lda <- align_components(Phi_lda, fit$Phi)
perm_stm <- align_components(Phi_stm, fit$Phi)
Theta_lda_a <- Theta_lda[, perm_lda, drop = FALSE]
Theta_stm_a <- Theta_stm[, perm_stm, drop = FALSE]

## ---------------------------------------------------------------------------
## 4. ALR-WLS covariate effects for LDA and STM
## ---------------------------------------------------------------------------
res_lda <- .alr_wls(Theta_lda_a, X_std, ref = K)
res_stm <- .alr_wls(Theta_stm_a, X_std, ref = K)

.tidy_effects <- function(res, method_tag) {
  P <- nrow(res$B); K1 <- ncol(res$B)
  data.frame(
    method   = method_tag,
    covariate = rep(rownames(res$B), times = K1),
    component = rep(colnames(res$B), each = P),
    estimate  = as.numeric(res$B),
    std.error = as.numeric(res$se),
    conf.low  = as.numeric(res$ci_lo),
    conf.high = as.numeric(res$ci_hi),
    stringsAsFactors = FALSE
  ) %>%
    mutate(z = estimate / pmax(std.error, .Machine$double.eps),
           p.value = 2 * stats::pnorm(-abs(z)),
           p.adj   = stats::p.adjust(p.value, method = "BH"))
}

eff_lda <- .tidy_effects(res_lda, "LDA+ALR")
eff_stm <- .tidy_effects(res_stm, "STM")

## GSCA-MM effects from the fit (using covariate_effects)
eff_gsca_obj <- covariate_effects(fit, ref = K, level = 0.95, adjust = "BH")
co_g <- eff_gsca_obj$coefficients
co_g <- co_g[co_g$covariate != "(Intercept)", ]
eff_gsca <- data.frame(
  method   = "GSCA-MM",
  covariate = co_g$covariate,
  component = co_g$component,
  estimate  = co_g$estimate,
  std.error = co_g$std.error,
  conf.low  = co_g$conf.low,
  conf.high = co_g$conf.high,
  z         = co_g$statistic,
  p.value   = co_g$p.value,
  p.adj     = co_g$p.adj,
  stringsAsFactors = FALSE
)

eff_all <- bind_rows(eff_gsca, eff_lda, eff_stm)
write.csv(eff_lda,  file.path(RESULT_DIR, "effects_lda.csv"),  row.names = FALSE)
write.csv(eff_stm,  file.path(RESULT_DIR, "effects_stm.csv"),  row.names = FALSE)
write.csv(eff_all,  file.path(RESULT_DIR, "effects_all.csv"),  row.names = FALSE)

## ---------------------------------------------------------------------------
## 5. Comparison plot: dot plot per (covariate, topic), paneled by covariate,
##    with three colors for the three methods.
## ---------------------------------------------------------------------------
.cov_label <- c(
  "year2020"                       = "Year = 2020 (vs 2019)",
  "year2021"                       = "Year = 2021 (vs 2019)",
  "partydemocratic"                = "Party = Democratic (vs Republican)",
  "candidate_typevice_president"   = "Candidate = VP (vs President)"
)
eff_all <- eff_all %>%
  mutate(
    covariate_label = factor(.cov_label[covariate],
                             levels = unname(.cov_label)),
    component_label = factor(component,
                             levels = paste0("comp", seq_len(K - 1L)),
                             labels = paste0("Topic ", seq_len(K - 1L))),
    method = factor(method, levels = c("GSCA-MM", "LDA+ALR", "STM"))
  )

p_cmp <- ggplot(eff_all,
                aes(x = estimate, y = component_label,
                    colour = method, group = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.55),
                 height = 0.25, linewidth = 0.5) +
  geom_point(size = 2, position = position_dodge(width = 0.55)) +
  facet_wrap(~ covariate_label, scales = "free", ncol = 2) +
  scale_colour_manual(values = METHOD_PALETTE, name = NULL) +
  labs(title = "Covariate effects on topic prevalence -- method comparison",
       subtitle = "Topics aligned to GSCA-MM via the Hungarian algorithm on Phi",
       x = "ALR coefficient (log-odds vs reference topic)",
       y = NULL) +
  theme_paper()

ggsave(file.path(FIGURE_DIR, "method_comparison_effects.png"),
       plot = p_cmp,
       width = 11, height = 1.5 + 0.45 * (K - 1) * 2,
       dpi = 150, limitsize = FALSE)
cat(sprintf("Method comparison -> %s\n",
            file.path(FIGURE_DIR, "method_comparison_effects.png")))

## ---------------------------------------------------------------------------
## 6. Significant-effect count per (method, covariate)
## ---------------------------------------------------------------------------
sig_count <- eff_all %>%
  mutate(direction = case_when(
    p.adj < 0.05 & estimate > 0 ~ "Positive",
    p.adj < 0.05 & estimate < 0 ~ "Negative",
    TRUE                         ~ "Non-significant"
  )) %>%
  group_by(method, covariate_label, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(direction = factor(direction,
                            levels = c("Negative", "Non-significant",
                                       "Positive")))

p_sig <- ggplot(sig_count,
                aes(x = method, y = n, fill = direction)) +
  geom_col() +
  facet_wrap(~ covariate_label, ncol = 2) +
  scale_fill_manual(values = c("Positive"        = "#1E40AF",
                                "Negative"        = "#B91C1C",
                                "Non-significant" = "grey70"),
                    name = NULL) +
  labs(title = "Significant covariate effects across methods",
       subtitle = sprintf("Counts out of K-1 = %d topic-level coefficients per panel",
                          K - 1L),
       x = NULL, y = "Number of effects") +
  theme_paper()

ggsave(file.path(FIGURE_DIR, "method_comparison_signif_count.png"),
       plot = p_sig,
       width = 10, height = 6, dpi = 150)
cat(sprintf("Sig-count plot -> %s\n",
            file.path(FIGURE_DIR, "method_comparison_signif_count.png")))
