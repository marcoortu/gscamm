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
##
## The design is configurable via environment variables so the same script
## can drive the incremental ablation runs (B/C/D/F) without code edits:
##   GSCAMM_FIT_MAX_ITER     gscamm EM cap            default 200
##   GSCAMM_BOOT_MAX_ITER    bootstrap EM cap         default 60
##   GSCAMM_NOISE_SCALE      bootstrap noise scale    default 1.0
##   GSCAMM_USE_POLISH       gscamm MAP polish on/off default 1 (on)
##   GSCAMM_SIGMA2_POLISH    sigma^2 for MAP polish   default 0.25
##   GSCAMM_USE_STM_RANDOM   include stm_random      default 1 (on)
##   GSCAMM_SIM_R            replicates per scenario default 100
##   GSCAMM_BOOT_B           bootstrap reps          default 200
## ---------------------------------------------------------------------------

.env_int   <- function(key, default) {
  v <- Sys.getenv(key, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) return(as.integer(default))
  as.integer(v)
}
.env_num   <- function(key, default) {
  v <- Sys.getenv(key, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) return(as.numeric(default))
  as.numeric(v)
}
.env_bool  <- function(key, default) {
  v <- Sys.getenv(key, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) return(as.logical(default))
  v <- tolower(v)
  v %in% c("1", "true", "t", "yes", "y", "on")
}

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
  noise_scale      = .env_num("GSCAMM_NOISE_SCALE", 1.0),
  bias_correct     = TRUE,
  ## fit controls
  fit_max_iter   = .env_int("GSCAMM_FIT_MAX_ITER", 200),
  fit_tol        = 1e-5,
  boot_max_iter  = .env_int("GSCAMM_BOOT_MAX_ITER", 60),
  boot_tol       = 1e-4,
  ## polish (MAP) controls
  use_polish     = .env_bool("GSCAMM_USE_POLISH", TRUE),
  sigma2_polish  = .env_num("GSCAMM_SIGMA2_POLISH", 0.25),
  ## parametric bootstrap on the plug-in fit: produces an additional
  ## coverage column on the gscamm row (coverage_B_boot_param). Set to 0
  ## to skip; default 200 matches the legacy noise-augmented bootstrap.
  param_boot_B          = .env_int("GSCAMM_PARAM_BOOT_B", 200),
  param_boot_max_iter   = .env_int("GSCAMM_PARAM_BOOT_MAX_ITER", 30),
  ## comparator switches: STM is fit with Random init by default (model-
  ## to-model fair comparison). The Spectral (anchor-words) variant is an
  ## opt-in counterfactual that quantifies the warm-start contribution
  ## but is NOT a model comparison -- it is an init comparison.
  use_stm_spectral = .env_bool("GSCAMM_USE_STM_SPECTRAL", FALSE)
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
       time = as.numeric(Sys.time() - t0, units = "secs"),
       theta_kind = "posterior",
       phi_init = "random_dirichlet")
}

## ---- STM comparator ------------------------------------------------------
##
## Two STM variants are provided:
##   fit_stm()           init.type = "Random"    (no warm start; canonical)
##   fit_stm_spectral()  init.type = "Spectral"  (anchor-words; opt-in counter-factual)
## The Spectral init is a strong data-driven warm start (Arora-Halpern-Mimno
## anchor-words algorithm). The ablation analysis (full_metrics_run{C,D,F})
## showed it contributes ~70-200% of STM's apparent recovery edge over the
## other models -- it is an init effect, not a model effect. To keep the
## comparison apples-to-apples between models we use the Random init by
## default; Spectral is reported only when GSCAMM_USE_STM_SPECTRAL=1 as an
## explicit appendix-style counter-factual.
.fit_stm_impl <- function(W, X, K, ref = K, max_iter = 75,
                          init_type = c("Spectral", "Random", "LDA"),
                          ...) {
  init_type <- match.arg(init_type)
  t0 <- Sys.time()
  W_int <- round(W); storage.mode(W_int) <- "integer"
  V <- ncol(W_int)
  vocab <- if (is.null(colnames(W_int))) sprintf("v%d", seq_len(V))
           else colnames(W_int)

  documents <- lapply(seq_len(nrow(W_int)), function(i) {
    nz <- which(W_int[i, ] > 0)
    if (!length(nz)) return(matrix(0L, nrow = 2L, ncol = 0L))
    rbind(as.integer(nz), as.integer(W_int[i, nz]))
  })

  proc <- stm::prepDocuments(documents, vocab,
                             lower.thresh = 0L, verbose = FALSE)
  kept_idx <- match(proc$vocab, vocab)

  X_df <- as.data.frame(X)
  prev_form <- as.formula(paste("~",
    paste(colnames(X_df), collapse = "+")))
  fit <- stm::stm(documents = proc$documents, vocab = proc$vocab,
                  K = K, prevalence = prev_form, data = X_df,
                  max.em.its = max_iter, init.type = init_type,
                  verbose = FALSE)
  Theta <- fit$theta
  Phi_kept <- exp(fit$beta$logbeta[[1]])
  Phi_kept <- Phi_kept / rowSums(Phi_kept)
  Phi <- matrix(0, K, V)
  Phi[, kept_idx] <- Phi_kept
  rs <- rowSums(Phi); rs[rs == 0] <- 1; Phi <- Phi / rs

  X_std <- scale(X); attr(X_std, "scaled:center") <- NULL
  attr(X_std, "scaled:scale") <- NULL
  res <- .alr_wls(Theta, X_std, ref = ref)
  list(Theta = Theta, Phi = Phi,
       B_alr = res$B, se_B = res$se,
       ci_lo = res$ci_lo, ci_hi = res$ci_hi,
       perplexity = gscamm::perplexity(W, Theta, Phi),
       time = as.numeric(Sys.time() - t0, units = "secs"),
       theta_kind = "posterior",
       phi_init = tolower(init_type))
}

fit_stm <- function(W, X, K, ref = K, max_iter = 75, ...)
  .fit_stm_impl(W, X, K, ref = ref, max_iter = max_iter,
                init_type = "Random", ...)

fit_stm_spectral <- function(W, X, K, ref = K, max_iter = 75, ...)
  .fit_stm_impl(W, X, K, ref = ref, max_iter = max_iter,
                init_type = "Spectral", ...)

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
##
## fit_gscamm_pi() returns the plug-in WLS estimates AND, optionally, the
## parametric-bootstrap CIs on the same fit (param_boot_B > 0). The
## parametric bootstrap resamples W^(b) ~ Mult(L_i, theta_hat_i' Phi_hat),
## warm-starts each refit at (Phi_hat, B_hat), and produces a basic CI
## with implicit bias correction. The plug-in estimates returned by this
## function are bit-identical to the no-bootstrap path: the bootstrap
## adds CI fields but never modifies the plug-in B_alr / se_B / ci_lo /
## ci_hi or any other downstream point estimate.
fit_gscamm_pi <- function(W, X, K, ref = K,
                          fit_max_iter = 80, fit_tol = 1e-5,
                          link = "logistic_normal",
                          gsca_space = "alr",
                          polish = "map",
                          sigma2_polish = 0.25,
                          param_boot_B = 0L,
                          param_boot_max_iter = 30L,
                          param_boot_seed = NULL, ...) {
  t0 <- Sys.time()
  ctl <- gscamm_control(max_iter = fit_max_iter, tol = fit_tol,
                        polish_max_iter = 50L)
  fit <- fit_gscamm(W, X, K = K, link = link,
                    gsca_space = gsca_space, gsca_ref = ref,
                    polish = polish, sigma2_polish = sigma2_polish,
                    control = ctl, seed = NULL)
  ## plug-in WLS
  eff <- covariate_effects(fit, ref = ref)
  P <- ncol(X)
  B  <- eff$B_alr[2:(P + 1L), , drop = FALSE]
  lo <- matrix(eff$coefficients$conf.low,  P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]
  hi <- matrix(eff$coefficients$conf.high, P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]
  se <- matrix(eff$coefficients$std.error, P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]

  ## ----- parametric bootstrap (additive: never modifies plug-in fields)
  lo_boot <- NULL; hi_boot <- NULL; t_boot <- NA_real_
  if (param_boot_B > 0L) {
    tb0 <- Sys.time()
    ctl_b <- gscamm_control(max_iter = param_boot_max_iter, tol = fit_tol,
                            polish_max_iter = 1L)
    eff_pb <- tryCatch(
      bootstrap_covariate_effects(fit, B = param_boot_B,
                                  method = "parametric",
                                  type = "basic",
                                  level = 0.95,
                                  ref = ref,
                                  control = ctl_b,
                                  seed = param_boot_seed,
                                  parallel = FALSE,
                                  adjust = "none"),
      error = function(e) {
        message("param-boot failed: ", conditionMessage(e)); NULL })
    if (!is.null(eff_pb)) {
      lo_boot <- matrix(eff_pb$coefficients$conf.low,
                        P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]
      hi_boot <- matrix(eff_pb$coefficients$conf.high,
                        P + 1L, K - 1L)[2:(P + 1L), , drop = FALSE]
    }
    t_boot <- as.numeric(Sys.time() - tb0, units = "secs")
  }

  ## R is the data-informed token-level posterior; Theta_map is the per-row
  ## MAP of eta combining the structural prior X*B with the data likelihood.
  ## Theta_map is the apples-to-apples counterpart of LDA gamma / STM theta.
  R_post <- fit$R %||% fit$Theta
  list(Theta = fit$Theta, R = R_post, Theta_map = fit$Theta_map,
       Phi = fit$Phi,
       B_alr = B, se_B = se, ci_lo = lo, ci_hi = hi,
       ci_lo_param_boot = lo_boot, ci_hi_param_boot = hi_boot,
       perplexity = gscamm::perplexity(W, R_post, fit$Phi),
       time = as.numeric(Sys.time() - t0, units = "secs"),
       time_param_boot = t_boot,
       theta_kind = "structural",
       phi_init = fit$init_phi,
       iters_converged = fit$convergence$iterations,
       converged = fit$convergence$converged,
       polish_n_converged = if (!is.null(fit$convergence$polish$n_converged))
         fit$convergence$polish$n_converged else NA_integer_,
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
                            polish = "map",
                            sigma2_polish = 0.25,
                            reuse_fit = NULL, ...) {
  t0 <- Sys.time()
  if (is.null(reuse_fit)) {
    ctl <- gscamm_control(max_iter = fit_max_iter, tol = fit_tol,
                          polish_max_iter = 50L)
    fit <- fit_gscamm(W, X, K = K, link = link,
                      gsca_space = gsca_space, gsca_ref = ref,
                      polish = polish, sigma2_polish = sigma2_polish,
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
  list(Theta = fit$Theta, R = R_post, Theta_map = fit$Theta_map,
       Phi = fit$Phi,
       B_alr = B, se_B = se, ci_lo = lo, ci_hi = hi,
       perplexity = gscamm::perplexity(W, R_post, fit$Phi),
       time = as.numeric(Sys.time() - t0, units = "secs"),
       theta_kind = "structural",
       phi_init = fit$init_phi,
       iters_converged = fit$convergence$iterations,
       converged = fit$convergence$converged,
       polish_n_converged = if (!is.null(fit$convergence$polish$n_converged))
         fit$convergence$polish$n_converged else NA_integer_)
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
  ## MAP polish (gscamm only). For LDA/STM the posterior Theta IS already
  ## the prior+data combined estimate, so we equate Theta_map_a = Theta_a
  ## to keep the rmse_theta_map column comparable across methods.
  Theta_map_a <- if (!is.null(fit_out$Theta_map))
    fit_out$Theta_map[, perm, drop = FALSE] else Theta_a

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
    Theta_map_a <- if (!is.null(fit_out$Theta_map))
      fit_out$Theta_map[, perm2, drop = FALSE] else Theta_a
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

  ## parametric-bootstrap CIs (only present for gscamm plug-in path)
  if (!is.null(fit_out$ci_lo_param_boot)) {
    lo_pb <- fit_out$ci_lo_param_boot[, fit_nr_order, drop = FALSE]
    hi_pb <- fit_out$ci_hi_param_boot[, fit_nr_order, drop = FALSE]
  } else {
    lo_pb <- NULL; hi_pb <- NULL
  }

  ## metrics across the three theta estimators:
  ##   rmse_theta      (structural, deterministic for gscamm; posterior for LDA/STM)
  ##   rmse_theta_R    (token-level responsibility posterior for gscamm; same as
  ##                    rmse_theta for LDA/STM where R = Theta posterior)
  ##   rmse_theta_map  (per-row MAP combining structural prior with the
  ##                    multinomial likelihood; equals rmse_theta for LDA/STM)
  ## The MAP column is the apples-to-apples target: it lines up with what
  ## LDA reports as gamma and STM reports as theta -- a covariate-aware
  ## prior combined with the data likelihood.
  rmse_theta     <- sqrt(mean((Theta_a     - sim$Theta)^2))
  rmse_theta_R   <- sqrt(mean((R_a         - sim$Theta)^2))
  rmse_theta_map <- sqrt(mean((Theta_map_a - sim$Theta)^2))
  rmse_phi       <- sqrt(mean((Phi_a       - sim$Phi)^2))
  rmse_B         <- sqrt(mean((B_a - sim$beta)^2))
  cov_B   <- if (!is.null(lo_a))
    mean(sim$beta >= lo_a & sim$beta <= hi_a) else NA_real_
  width_B <- if (!is.null(lo_a)) mean(hi_a - lo_a) else NA_real_

  ## parametric-bootstrap coverage / width: populated only for the gscamm
  ## plug-in path; NA for gscamm_boot, lda, stm, stm_spectral.
  cov_B_boot_param   <- if (!is.null(lo_pb))
    mean(sim$beta >= lo_pb & sim$beta <= hi_pb) else NA_real_
  width_B_boot_param <- if (!is.null(lo_pb)) mean(hi_pb - lo_pb) else NA_real_

  list(method = method_tag,
       rmse_theta = rmse_theta,
       rmse_theta_R = rmse_theta_R,
       rmse_theta_map = rmse_theta_map,
       rmse_phi = rmse_phi,
       rmse_B = rmse_B, coverage_B = cov_B, width_B = width_B,
       coverage_B_boot_param = cov_B_boot_param,
       width_B_boot_param    = width_B_boot_param,
       perplexity = fit_out$perplexity, time = fit_out$time,
       time_param_boot = fit_out$time_param_boot %||% NA_real_,
       theta_kind = fit_out$theta_kind %||% NA_character_,
       phi_init   = fit_out$phi_init   %||% NA_character_,
       iters_converged   = fit_out$iters_converged   %||% NA_integer_,
       converged         = fit_out$converged         %||% NA,
       polish_n_converged = fit_out$polish_n_converged %||% NA_integer_)
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
  polish_arg <- if (isTRUE(design$use_polish)) "map" else "none"

  out <- list()
  ## gscamm plug-in (with optional parametric bootstrap CIs as a column,
  ## NOT as a separate row -- design$param_boot_B controls B; 0 disables)
  m1 <- fit_gscamm_pi(sim$W, sim$X, K, ref = K,
                      fit_max_iter = design$fit_max_iter,
                      fit_tol = design$fit_tol,
                      link = link_for_gscamm,
                      polish = polish_arg,
                      sigma2_polish = design$sigma2_polish,
                      param_boot_B = design$param_boot_B %||% 0L,
                      param_boot_max_iter = design$param_boot_max_iter %||% 30L)
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
                        polish = polish_arg,
                        sigma2_polish = design$sigma2_polish,
                        reuse_fit = m1$fit)
  out$gscamm_boot <- align_and_metrics(m2, sim, ref = K,
                                       method_tag = "gscamm_boot")

  ## LDA + ALR
  m3 <- tryCatch(fit_lda_alr(sim$W, sim$X, K, ref = K),
                 error = function(e) NULL)
  if (!is.null(m3))
    out$lda <- align_and_metrics(m3, sim, ref = K, method_tag = "lda")

  ## STM with Random init (canonical model-to-model comparison; native
  ## Spectral init is a strong anchor-words warm start that adds 35-200%
  ## init-driven advantage on theta recovery, so it would not be a fair
  ## model comparison -- see ablation runC/D/F vs stm_random in the paper).
  m4 <- tryCatch(fit_stm(sim$W, sim$X, K, ref = K),
                 error = function(e) { message("STM error: ", conditionMessage(e)); NULL })
  if (!is.null(m4))
    out$stm <- align_and_metrics(m4, sim, ref = K, method_tag = "stm")

  ## STM with Spectral (anchor-words) init: opt-in counter-factual that
  ## quantifies the warm-start contribution to STM's apparent recovery edge.
  ## Reported as an explicit init-comparison row, not as a model alternative.
  if (isTRUE(design$use_stm_spectral)) {
    m5 <- tryCatch(fit_stm_spectral(sim$W, sim$X, K, ref = K),
                   error = function(e) {
                     message("STM(spectral) error: ", conditionMessage(e)); NULL })
    if (!is.null(m5))
      out$stm_spectral <- align_and_metrics(m5, sim, ref = K,
                                            method_tag = "stm_spectral")
  }

  do.call(rbind, lapply(out, function(o) data.frame(o, scenario = scenario,
                                                    stringsAsFactors = FALSE)))
}
