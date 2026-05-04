## ---------------------------------------------------------------------------
## 00_setup.R -- common settings, paths, helper functions used across the
## simulation scripts (01_run_pilot.R, 02_run_full.R, 03_make_table1.R,
## 04_make_figures.R).
##
## Run from the project root.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(gscamm)
  library(topicmodels)   ## LDA
  library(stm)           ## STM
  library(slam)          ## simple_triplet_matrix used by topicmodels
})

## null-coalescing helper (gscamm's `%||%` is not exported)
`%||%` <- function(a, b) if (is.null(a)) b else a

## resolve paths relative to this file's location
.this_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).*", args, perl = TRUE))
  if (length(m)) normalizePath(m[1])
  else if (!is.null(sys.frames())) sys.frames()[[1]]$ofile
  else NA_character_
}
HERE <- tryCatch(dirname(.this_file()), error = function(e) getwd())
if (!dir.exists(HERE)) HERE <- file.path(getwd(),
  "replication_package/simulations/code")

REPL_ROOT  <- normalizePath(file.path(HERE, ".."))   ## .../simulations
RESULT_DIR <- file.path(REPL_ROOT, "results")
FIGURE_DIR <- file.path(REPL_ROOT, "figures")
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)
if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

## ---------------------------------------------------------------------------
## Default Monte Carlo design (paper: N=1000, V=500, K=10, P=8, R=100).
## Override in the calling script if you want a smaller pilot.
## ---------------------------------------------------------------------------
DEFAULT_DESIGN <- list(
  N = 1000, V = 500, K = 10, P = 8,
  R = 100,
  scenarios = c("baseline", "high_covariate", "high_sparsity"),
  links = c(baseline = "logistic_normal",
            high_covariate = "dirichlet",
            high_sparsity  = "zero_inflated"),
  ## bootstrap settings for the gscamm_boot method
  boot_B           = 200,
  boot_type        = "basic",
  noise_augment    = TRUE,
  noise_scale      = 2.0,
  bias_correct     = TRUE,
  ## fit controls
  fit_max_iter   = 80,
  fit_tol        = 1e-5,
  boot_max_iter  = 30,
  boot_tol       = 1e-4
)

## ---------------------------------------------------------------------------
## Helper: shared fit_method dispatch returning a uniform list.
## Returns:
##   $Theta   N x K mixture proportions
##   $Phi     K x V component-category distributions
##   $B_alr   P x (K-1) ALR coefficient matrix (slopes, no intercept row)
##   $se_B    P x (K-1) standard errors aligned with B_alr (or NULL)
##   $ci_lo, $ci_hi  P x (K-1) CI bounds (or NULL if SE is NULL)
##   $perplexity scalar
##   $time    elapsed seconds
## ---------------------------------------------------------------------------

## ---- LDA + ALR-WLS comparator --------------------------------------------
fit_lda_alr <- function(W, X, K, ref = K, alpha_lda = 0.1, ...) {
  t0 <- Sys.time()
  ## topicmodels expects integer counts; convert to triplet for memory
  W_int <- round(W)
  storage.mode(W_int) <- "integer"
  dtm <- slam::as.simple_triplet_matrix(W_int)

  ## fit LDA via Variational EM (Blei et al., 2003)
  ld <- topicmodels::LDA(dtm, k = K,
                         method = "VEM",
                         control = list(alpha = alpha_lda, seed = 1))
  Phi   <- exp(ld@beta)                          ## K x V
  Phi   <- Phi / rowSums(Phi)
  Theta <- ld@gamma                              ## N x K (posterior means)

  ## ALR-WLS on standardized X
  X_std <- scale(X)
  attr(X_std, "scaled:center") <- NULL
  attr(X_std, "scaled:scale")  <- NULL
  res <- .alr_wls(Theta, X_std, ref = ref)

  list(Theta = Theta, Phi = Phi,
       B_alr = res$B, se_B = res$se,
       ci_lo = res$ci_lo, ci_hi = res$ci_hi,
       perplexity = gscamm::perplexity(W, Theta, Phi),
       time = as.numeric(Sys.time() - t0, units = "secs"))
}

## ---- STM comparator ------------------------------------------------------
fit_stm <- function(W, X, K, ref = K, max_iter = 75, ...) {
  t0 <- Sys.time()
  W_int <- round(W); storage.mode(W_int) <- "integer"
  V <- ncol(W_int)
  vocab <- if (is.null(colnames(W_int))) sprintf("v%d", seq_len(V))
           else colnames(W_int)

  ## convert each row of W into stm's (term, count) integer matrix
  documents <- lapply(seq_len(nrow(W_int)), function(i) {
    nz <- which(W_int[i, ] > 0)
    if (!length(nz)) return(matrix(0L, nrow = 2L, ncol = 0L))
    rbind(as.integer(nz), as.integer(W_int[i, nz]))
  })

  ## stm requires every vocab term to appear in at least one document.
  ## prepDocuments drops zero-count terms and remaps indices; we keep
  ## track of the surviving indices to re-expand Phi back to length V.
  proc <- stm::prepDocuments(documents, vocab,
                             lower.thresh = 0L, verbose = FALSE)
  kept_idx <- match(proc$vocab, vocab)               ## original V positions

  X_df <- as.data.frame(X)
  prev_form <- as.formula(paste("~",
    paste(colnames(X_df), collapse = "+")))
  fit <- stm::stm(documents = proc$documents, vocab = proc$vocab,
                  K = K, prevalence = prev_form, data = X_df,
                  max.em.its = max_iter, init.type = "Spectral",
                  verbose = FALSE)
  Theta <- fit$theta                                  ## N x K
  Phi_kept <- exp(fit$beta$logbeta[[1]])              ## K x length(kept_idx)
  Phi_kept <- Phi_kept / rowSums(Phi_kept)
  ## re-expand to full V (dropped terms get 0 weight)
  Phi <- matrix(0, K, V)
  Phi[, kept_idx] <- Phi_kept
  ## guard against rows that are now all zero (shouldn't happen unless
  ## all top words for a topic were dropped)
  rs <- rowSums(Phi); rs[rs == 0] <- 1; Phi <- Phi / rs

  X_std <- scale(X); attr(X_std, "scaled:center") <- NULL
  attr(X_std, "scaled:scale") <- NULL
  res <- .alr_wls(Theta, X_std, ref = ref)
  list(Theta = Theta, Phi = Phi,
       B_alr = res$B, se_B = res$se,
       ci_lo = res$ci_lo, ci_hi = res$ci_hi,
       perplexity = gscamm::perplexity(W, Theta, Phi),
       time = as.numeric(Sys.time() - t0, units = "secs"))
}

## ---- shared ALR-WLS used by LDA and STM comparators ----------------------
.alr_wls <- function(Theta, X_std, ref = ncol(Theta), level = 0.95) {
  Theta <- pmax(Theta, .Machine$double.eps)
  N <- nrow(Theta); K <- ncol(Theta); P <- ncol(X_std)
  D <- cbind(`(Intercept)` = 1, X_std)
  non_ref <- setdiff(seq_len(K), ref)
  K1 <- length(non_ref)
  B <- matrix(NA_real_, P, K1)
  se <- matrix(NA_real_, P, K1)
  z <- stats::qnorm(1 - (1 - level) / 2)
  for (j in seq_along(non_ref)) {
    k <- non_ref[j]
    y <- log(Theta[, k] / Theta[, ref])
    w <- Theta[, k] * Theta[, ref]
    f <- stats::lm.wfit(D, y, w)
    bhat <- f$coefficients
    bhat[is.na(bhat)] <- 0
    e <- as.numeric(f$residuals)
    XtWX <- crossprod(D, w * D)
    XtWX_inv <- tryCatch(solve(XtWX), error = function(...) NULL)
    if (is.null(XtWX_inv)) next
    meat <- crossprod(D, (w^2 * e^2) * D)
    Vmat <- XtWX_inv %*% meat %*% XtWX_inv
    B[, j]  <- bhat[-1L]                           ## drop intercept row
    se[, j] <- sqrt(pmax(diag(Vmat)[-1L], 0))
  }
  list(B = B, se = se,
       ci_lo = B - z * se, ci_hi = B + z * se)
}

## ---- gscamm comparators (plug-in and bootstrap variants) -----------------
fit_gscamm_pi <- function(W, X, K, ref = K,
                          fit_max_iter = 80, fit_tol = 1e-5,
                          link = "logistic_normal",
                          gsca_space = "alr", ...) {
  t0 <- Sys.time()
  ctl <- gscamm_control(max_iter = fit_max_iter, tol = fit_tol)
  fit <- fit_gscamm(W, X, K = K, link = link,
                    gsca_space = gsca_space, gsca_ref = ref,
                    control = ctl, seed = NULL)
  ## plug-in WLS
  eff <- covariate_effects(fit, ref = ref)
  P <- ncol(X)
  B  <- eff$B_alr[2:(P + 1L), , drop = FALSE]
  lo <- matrix(eff$coefficients$conf.low,  P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]
  hi <- matrix(eff$coefficients$conf.high, P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]
  se <- matrix(eff$coefficients$std.error, P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]

  ## R is the data-informed posterior (analog of LDA's gamma and STM's
  ## theta); the deterministic fit$Theta is the covariate-only prior.
  ## We expose R so align_and_metrics can compare it to sim$Theta on the
  ## same footing as the LDA/STM comparators, and compute perplexity from
  ## R for consistency with how those packages report it.
  R_post <- fit$R %||% fit$Theta
  list(Theta = fit$Theta, R = R_post, Phi = fit$Phi,
       B_alr = B, se_B = se, ci_lo = lo, ci_hi = hi,
       perplexity = gscamm::perplexity(W, R_post, fit$Phi),
       time = as.numeric(Sys.time() - t0, units = "secs"),
       fit = fit)                                   ## kept for reuse by boot
}

fit_gscamm_boot <- function(W, X, K, ref = K,
                            fit_max_iter = 80, fit_tol = 1e-5,
                            boot_max_iter = 30, boot_tol = 1e-4,
                            boot_B = 200, boot_type = "basic",
                            noise_augment = TRUE, noise_scale = 2.0,
                            bias_correct = TRUE,
                            link = "logistic_normal",
                            gsca_space = "alr",
                            reuse_fit = NULL, ...) {
  t0 <- Sys.time()
  if (is.null(reuse_fit)) {
    ctl <- gscamm_control(max_iter = fit_max_iter, tol = fit_tol)
    fit <- fit_gscamm(W, X, K = K, link = link,
                      gsca_space = gsca_space, gsca_ref = ref,
                      control = ctl, seed = NULL)
  } else {
    fit <- reuse_fit
  }
  ctl_b <- gscamm_control(max_iter = boot_max_iter, tol = boot_tol)
  eff <- bootstrap_covariate_effects(
    fit, B = boot_B, ref = ref, type = boot_type,
    noise_augment = noise_augment, noise_scale = noise_scale,
    bias_correct = bias_correct,
    control = ctl_b, seed = NULL)
  P <- ncol(X)
  B  <- eff$B_alr[2:(P + 1L), , drop = FALSE]
  lo <- matrix(eff$coefficients$conf.low,  P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]
  hi <- matrix(eff$coefficients$conf.high, P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]
  se <- matrix(eff$coefficients$std.error, P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]
  R_post <- fit$R %||% fit$Theta
  list(Theta = fit$Theta, R = R_post, Phi = fit$Phi,
       B_alr = B, se_B = se, ci_lo = lo, ci_hi = hi,
       perplexity = gscamm::perplexity(W, R_post, fit$Phi),
       time = as.numeric(Sys.time() - t0, units = "secs"))
}

## ---------------------------------------------------------------------------
## Alignment + metric computation. Aligns a fitted (Theta_hat, Phi_hat,
## B_alr) to the truth via Hungarian on Phi, then computes RMSE and
## coverage. Crucially, the alignment must respect the ALR reference
## convention: the truth uses comp K as reference, so after alignment
## comp K of the fit must match the true reference.
## ---------------------------------------------------------------------------
align_and_metrics <- function(fit_out, sim, ref = ncol(sim$Phi),
                              method_tag = "method") {
  K <- nrow(sim$Phi); P <- ncol(sim$X)

  ## permutation: align fit$Phi rows to sim$Phi rows
  perm <- gscamm::align_components(fit_out$Phi, sim$Phi)
  Phi_a   <- fit_out$Phi[perm, , drop = FALSE]
  Theta_a <- fit_out$Theta[, perm, drop = FALSE]
  ## data-informed posterior (R = posterior responsibilities for gscamm,
  ## fit_out$Theta itself for LDA/STM where Theta is already data-informed).
  R_a <- if (!is.null(fit_out$R))
    fit_out$R[, perm, drop = FALSE] else Theta_a

  ## permute B_alr columns: B_alr columns correspond to non-reference
  ## components in the FIT's order (1..K-1). The fit's reference is comp
  ## K; after the permutation comp K_aligned corresponds to row perm[K]
  ## of the original fit. We need the COLUMN permutation that takes
  ## fit's non-reference-column ordering into the truth's non-reference
  ## ordering (also 1..K-1 with comp K as ref).
  ##
  ## fit's non-ref columns are at positions 1..K-1 in fit-order.
  ## After perm: fit_aligned[, k] = fit[, perm[k]].
  ## Non-ref columns of fit_aligned correspond to original positions
  ## perm[1], ..., perm[K-1]. We need to find which of fit's original
  ## non-ref columns map to which aligned non-ref column.
  ##
  ## Simpler: compute B_alr_aligned by re-running ALR-WLS on Theta_a if
  ## B_alr is not directly transformable.

  ## Re-run ALR-WLS on the aligned Theta to recover the aligned B_alr
  ## using the same standardized X. This avoids ambiguity.
  X_std <- scale(sim$X)
  attr(X_std, "scaled:center") <- NULL
  attr(X_std, "scaled:scale")  <- NULL

  ## If the method already supplied CIs, we want to permute them rather
  ## than recompute (because plug-in vs bootstrap CIs differ). The CIs
  ## are tied to fit's non-reference-column ordering. To re-align them,
  ## we need the column permutation.
  ## fit's non-ref columns are 1..K-1 (assuming ref = K).
  ## After permutation perm, the aligned columns are perm.
  ## Aligned non-ref columns are setdiff(seq_len(K), which(perm == K)).
  ## i.e., the position of perm == K in the aligned vector is the new
  ## reference position.
  new_ref_pos <- which(perm == ref)
  if (new_ref_pos != ref) {
    ## need to ALSO swap the new_ref_pos column with the K-th to put
    ## reference back at position K. This is a second permutation.
    swap <- seq_len(K)
    swap[new_ref_pos] <- ref
    swap[ref] <- new_ref_pos
    perm2 <- perm[swap]
    Phi_a   <- fit_out$Phi[perm2, , drop = FALSE]
    Theta_a <- fit_out$Theta[, perm2, drop = FALSE]
    R_a <- if (!is.null(fit_out$R))
      fit_out$R[, perm2, drop = FALSE] else Theta_a
  }
  ## At this point comp K of (Theta_a, Phi_a) is aligned with the true ref.

  ## Permute the supplied CIs accordingly. fit_out$B_alr/se_B/ci_lo/ci_hi
  ## are P x (K-1), columns in fit-order non-ref (1..K-1).
  ## We need to map each true non-ref index to the aligned non-ref index.
  fit_nr_order <- if (exists("perm2")) perm2[-K] else perm[-K]
  ## fit_nr_order[k] tells which original column-of-fit ends up at
  ## aligned position k (for k = 1..K-1, all non-ref).
  ## So aligned B[, k] = fit_out$B_alr[, idx(fit_nr_order[k])]
  ## where idx() maps original-fit-column j to its column in the
  ## (P x (K-1)) supplied B_alr. Original non-ref columns of fit are
  ## also 1..K-1 (since fit's ref is K). So idx is identity for j != K.
  ## fit_nr_order entries are guaranteed != K because we swapped K to
  ## the last position. Hence we can index directly.
  if (!is.null(fit_out$B_alr)) {
    B_a   <- fit_out$B_alr[, fit_nr_order, drop = FALSE]
    if (!is.null(fit_out$ci_lo)) {
      lo_a <- fit_out$ci_lo[, fit_nr_order, drop = FALSE]
      hi_a <- fit_out$ci_hi[, fit_nr_order, drop = FALSE]
    } else { lo_a <- NULL; hi_a <- NULL }
    if (!is.null(fit_out$se_B)) {
      se_a <- fit_out$se_B[, fit_nr_order, drop = FALSE]
    } else { se_a <- NULL }
  } else {
    res <- .alr_wls(Theta_a, X_std, ref = K)
    B_a <- res$B; lo_a <- res$ci_lo; hi_a <- res$ci_hi; se_a <- res$se
  }

  ## metrics: rmse_theta uses the structural estimator Theta_a (= deterministic
  ## inverse-ALR(X*B_hat) for gscamm; data-informed gamma/theta for LDA/STM).
  ## Under GSCA-MM's deterministic specification, Theta IS the Bayesian
  ## estimator of theta. The responsibilities R are the posterior over
  ## token-level assignments, NOT the posterior of theta, so they over-peak
  ## with informative data and underestimate the continuous fractional
  ## weights of the truth. We additionally report rmse_theta_R as a
  ## diagnostic for the data-informed responsibility-based estimator.
  rmse_theta   <- sqrt(mean((Theta_a - sim$Theta)^2))
  rmse_theta_R <- sqrt(mean((R_a     - sim$Theta)^2))
  rmse_phi     <- sqrt(mean((Phi_a   - sim$Phi)^2))
  rmse_B       <- sqrt(mean((B_a - sim$beta)^2))
  cov_B <- if (!is.null(lo_a))
    mean(sim$beta >= lo_a & sim$beta <= hi_a) else NA_real_
  width_B <- if (!is.null(lo_a)) mean(hi_a - lo_a) else NA_real_

  list(method = method_tag,
       rmse_theta = rmse_theta, rmse_theta_R = rmse_theta_R,
       rmse_phi = rmse_phi,
       rmse_B = rmse_B, coverage_B = cov_B, width_B = width_B,
       perplexity = fit_out$perplexity, time = fit_out$time)
}

## ---------------------------------------------------------------------------
## Run all four methods on a single (sim, scenario) and return a long
## data frame of metrics.
## ---------------------------------------------------------------------------
run_one_replicate <- function(sim, scenario, design = DEFAULT_DESIGN,
                              link_for_gscamm = NULL,
                              verbose = FALSE) {
  if (is.null(link_for_gscamm))
    link_for_gscamm <- design$links[[scenario]]
  K <- design$K

  out <- list()
  ## gscamm plug-in
  m1 <- fit_gscamm_pi(sim$W, sim$X, K, ref = K,
                      fit_max_iter = design$fit_max_iter,
                      fit_tol = design$fit_tol,
                      link = link_for_gscamm)
  out$gscamm <- align_and_metrics(m1, sim, ref = K, method_tag = "gscamm")

  ## gscamm bootstrap (reuse the same fit to halve the cost)
  m2 <- fit_gscamm_boot(sim$W, sim$X, K, ref = K,
                        boot_max_iter = design$boot_max_iter,
                        boot_tol = design$boot_tol,
                        boot_B = design$boot_B,
                        boot_type = design$boot_type,
                        noise_augment = design$noise_augment,
                        noise_scale = design$noise_scale,
                        bias_correct = design$bias_correct,
                        link = link_for_gscamm,
                        reuse_fit = m1$fit)
  out$gscamm_boot <- align_and_metrics(m2, sim, ref = K,
                                       method_tag = "gscamm_boot")

  ## LDA + ALR
  m3 <- tryCatch(fit_lda_alr(sim$W, sim$X, K, ref = K),
                 error = function(e) NULL)
  if (!is.null(m3))
    out$lda <- align_and_metrics(m3, sim, ref = K, method_tag = "lda")

  ## STM
  m4 <- tryCatch(fit_stm(sim$W, sim$X, K, ref = K),
                 error = function(e) { message("STM error: ", conditionMessage(e)); NULL })
  if (!is.null(m4))
    out$stm <- align_and_metrics(m4, sim, ref = K, method_tag = "stm")

  do.call(rbind, lapply(out, function(o) data.frame(o, scenario = scenario,
                                                    stringsAsFactors = FALSE)))
}
