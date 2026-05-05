## ---------------------------------------------------------------------------
## Bootstrap confidence intervals for the covariate-effect estimator B.
##
## Two resampling strategies are supported:
##
## (a) "parametric" -- under the deterministic model
##       theta_i = inverse-ALR(X_i' B_{-K})
##       y_i ~ Mult(L_i, theta_i' Phi)
##     the only stochastic source is multinomial token sampling. The
##     parametric bootstrap fixes (B_hat, Phi_hat) at the fitted values,
##     computes q_hat_i = theta_hat_i' Phi_hat, and resamples
##         W^(b)_i ~ Mult(L_i, q_hat_i)
##     for b = 1..B. Each (W^(b), X) is refit with warm-start at
##     (Phi_hat, B_hat); after Hungarian alignment to the original Phi
##     the basic-bootstrap CI is
##         [2 B_hat - q_{1-alpha/2},  2 B_hat - q_{alpha/2}]
##     which is implicitly bias-corrected.
##
## (b) "noise_augmented" -- the legacy non-parametric row bootstrap with
##     Gaussian noise injected into the ALR-space scores after each refit.
##     Designed for the (now superseded) generative formulation that
##     included a latent epsilon_ik term; retained for paper reproduction
##     and for inference under non-deterministic links.
##
## Reference: Ortu and Frigau (2026), Remark 2 and Section 4.2. The plug-in
## sandwich variance treats the first-stage estimate of Theta as if it
## were the truth and therefore undercovers when N is small relative to K
## or observation lengths are short -- the parametric bootstrap recovers
## near-nominal coverage for the path coefficients.
## ---------------------------------------------------------------------------

#' Bootstrap confidence intervals for covariate effects
#'
#' Implements the non-parametric row bootstrap for the path coefficients
#' estimated by GSCA-MM, propagating the first-stage uncertainty in
#' \eqn{\hat{\boldsymbol{\Theta}}} through to the second-stage ALR-WLS
#' regression. See Remark 2 of Ortu and Frigau (2026).
#'
#' @param fit a fitted \code{gscamm} object.
#' @param B number of bootstrap replications (default 200).
#' @param level confidence level for percentile intervals (default 0.95).
#' @param ref reference component for the ALR transform (default the last
#'   component of \code{fit}).
#' @param adjust multiple-testing adjustment for bootstrap p-values
#'   (default \code{"BH"}; see \code{\link[stats]{p.adjust}}).
#' @param control optional override of the \code{\link{gscamm_control}}
#'   used for the bootstrap re-fits. Defaults to \code{fit$control}; in
#'   practice it is sometimes useful to lower \code{max_iter} or relax
#'   \code{tol} to keep the bootstrap tractable.
#' @param seed integer for reproducibility.
#' @param verbose logical, print progress every 10 replicates.
#' @param parallel logical (default \code{FALSE}); if \code{TRUE}, uses
#'   \code{parallel::mclapply} on Unix or \code{parallel::parLapply} on
#'   Windows.
#' @param n_cores number of cores for parallel execution; defaults to
#'   half of \code{parallel::detectCores()}.
#' @param progress_callback optional function called as
#'   \code{progress_callback(b, B)} after each replicate completes; useful
#'   for embedding progress reporting in higher-level scripts.
#' @param type confidence-interval type. One of:
#'   \describe{
#'     \item{\code{"percentile"}}{(default) empirical quantiles of the
#'       bootstrap distribution.}
#'     \item{\code{"basic"}}{the basic / pivotal bootstrap CI,
#'       \eqn{[2\hat\theta - q_{1-\alpha/2}, 2\hat\theta - q_{\alpha/2}]}.}
#'     \item{\code{"bca"}}{bias-corrected and accelerated bootstrap with
#'       the acceleration constant computed via leave-one-out jackknife
#'       on the second-stage ALR-WLS estimator (the first-stage fit is
#'       not re-run for the jackknife, which keeps the cost negligible).}
#'   }
#' @param bias_correct logical (default \code{FALSE}). If \code{TRUE},
#'   shifts both the point estimate and the confidence interval by the
#'   bootstrap bias \eqn{\hat\theta - \bar{\theta}^*}. This is a
#'   first-order bias correction; combine with \code{type = "percentile"}
#'   for a clean bias-corrected percentile CI.
#' @param use response basis passed through to
#'   \code{\link{covariate_effects}}: \code{"theta"} (default) or
#'   \code{"responsibilities"}.
#' @param method resampling strategy. \code{"parametric"} (default) draws
#'   \eqn{W_i^{(b)} \sim \mathrm{Mult}(L_i, \hat{q}_i)} with
#'   \eqn{\hat{q}_i = \hat{\theta}_i^\top \hat{\Phi}} and refits with
#'   warm-start at \eqn{(\hat{\Phi}, \hat{B})}; under the deterministic
#'   GSCA-MM model this is the well-specified bootstrap because the only
#'   stochastic source is multinomial token sampling. \code{"noise_augmented"}
#'   is the legacy non-parametric row bootstrap with Gaussian noise injection
#'   in ALR space, provided for paper reproduction and for the (superseded)
#'   generative model that included a latent \eqn{\varepsilon_{ik}} term.
#' @param noise_augment legacy logical flag (default \code{FALSE}). When
#'   \code{TRUE} and \code{method} is not explicitly set, \code{method} is
#'   silently switched to \code{"noise_augmented"} for backward compatibility
#'   with code written against earlier versions of the package. When
#'   \code{method = "parametric"} this argument is ignored.
#'
#' @return an object of class \code{gscamm_effects} with the same shape as
#'   the output of \code{\link{covariate_effects}}, but with confidence
#'   intervals (\code{conf.low}/\code{conf.high}, \code{or.low}/\code{or.high})
#'   and \code{std.error} obtained from the bootstrap distribution
#'   (sd across replicates), and \code{p.value} a percentile-based two-sided
#'   p-value computed as \eqn{2 \min(F^*(0), 1 - F^*(0))} from the bootstrap
#'   ECDF of each coefficient. The list also contains:
#'   \item{B_draws}{\eqn{B \times (P+1) \times (K-1)} array of bootstrap
#'     ALR coefficients.}
#'   \item{B_draws_failed}{integer count of bootstrap replications that
#'     failed and were discarded.}
#' @export
bootstrap_covariate_effects <- function(fit,
                                        B = 200,
                                        level = 0.95,
                                        ref = fit$K,
                                        adjust = c("BH", "holm", "hochberg",
                                                   "hommel", "bonferroni",
                                                   "BY", "fdr", "none"),
                                        control = NULL,
                                        seed = NULL,
                                        verbose = FALSE,
                                        parallel = FALSE,
                                        n_cores = NULL,
                                        progress_callback = NULL,
                                        method = c("parametric",
                                                   "noise_augmented"),
                                        type = c("percentile", "basic", "bca"),
                                        bias_correct = FALSE,
                                        use = c("theta", "responsibilities"),
                                        noise_augment = FALSE,
                                        noise_scale = 1) {
  ## Back-compat: callers that pass noise_augment=TRUE without an explicit
  ## method get the legacy noise_augmented branch automatically.
  method_explicit <- !missing(method)
  type_explicit   <- !missing(type)
  method <- match.arg(method)
  type   <- match.arg(type)
  use    <- match.arg(use)
  if (!method_explicit && isTRUE(noise_augment)) method <- "noise_augmented"
  ## Parametric defaults to basic CIs (implicit bias correction); user can
  ## still pick percentile or BCa via the `type` argument.
  if (method == "parametric" && !type_explicit) type <- "basic"
  if (!inherits(fit, "gscamm"))
    stop("`fit` must be a gscamm object.")
  adjust <- match.arg(adjust)
  if (B < 10) stop("B must be >= 10 for percentile intervals.")
  if (is.null(control)) control <- fit$control
  if (is.null(fit$X))
    stop("`fit` does not contain the original covariate matrix; refit with the current package version.")

  ## prepare the original data
  W <- fit$W
  X <- fit$X
  N <- nrow(W)
  K <- fit$K
  P <- fit$P
  link <- fit$link
  ref_orig <- ref

  ## reference for shape: estimates from the original fit
  eff_orig <- covariate_effects(fit, ref = ref, level = level,
                                adjust = adjust, use = use)
  cov_names <- rownames(eff_orig$B_alr)
  comp_names <- colnames(eff_orig$B_alr)

  ## per-component residual sd in ALR space, estimated from the original
  ## responsibilities R via OLS of log(R[k]/R[ref]) on X. Used by the
  ## noise_augmented branch to recover the latent noise that a deterministic
  ## link discards.
  noise_sigma <- numeric(K - 1L)
  if (method == "noise_augmented") {
    Rorig <- pmax(fit$R %||% fit$Theta, .Machine$double.eps)
    Xs    <- fit$X_std
    Dols  <- cbind(1, Xs)
    nr_idx <- setdiff(seq_len(K), ref)
    for (j in seq_along(nr_idx)) {
      k <- nr_idx[j]
      y_alr <- log(Rorig[, k] / Rorig[, ref])
      f <- stats::lm.fit(Dols, y_alr)
      r <- as.numeric(f$residuals)
      df <- max(length(r) - ncol(Dols), 1L)
      noise_sigma[j] <- sqrt(sum(r^2) / df)
    }
    noise_sigma <- noise_sigma * noise_scale
  }

  ## --- precompute objects needed by the parametric branch -----------------
  ## Fitted multinomial probabilities q_hat_i = theta_hat_i' Phi_hat and the
  ## per-row token totals L_i. The parametric refits warm-start at
  ## (Phi_hat, B_hat) so they typically converge in <= 10 EM iterations;
  ## we cap max_iter accordingly to keep the bootstrap tractable.
  Theta_hat <- pmax(fit$Theta, .Machine$double.eps)
  Phi_hat   <- pmax(fit$Phi,   .Machine$double.eps)
  B_hat     <- fit$B
  q_hat     <- Theta_hat %*% Phi_hat                 ## N x V
  q_hat     <- pmax(q_hat, 0)
  q_hat     <- q_hat / pmax(rowSums(q_hat), .Machine$double.xmin)
  L_obs     <- as.integer(round(rowSums(W)))

  ## Refit control: NO warm-start at Phi or B.
  ##
  ## Validation pilot (2026-05-05) showed that warm-starting Phi at the
  ## fitted value collapses the bootstrap variance: with Phi pinned and W^(b)
  ## close to W, the M-step barely moves Phi, so B^(b) lands within ~3% of
  ## B_hat and the resulting CIs are 5-20x too narrow (coverage ~0.17,
  ## target 0.95). Each parametric refit must therefore go through the
  ## full EM with k-means++ initialization on W^(b) to reflect the true
  ## first-stage uncertainty. The cost overhead is modest because gscamm
  ## EM converges in O(50) iterations even from cold init when W^(b) is
  ## drawn from the fitted multinomial.
  control_boot <- control
  control_boot$init_Phi <- NULL
  control_boot$init_B   <- NULL
  control_boot$max_iter <- min(control$max_iter, 80L)

  ## one-shot worker
  one_boot <- function(b) {
    ## per-replicate seed (reproducible across cores when base seed given)
    if (!is.null(seed)) set.seed(seed + 7919L * b)

    if (method == "parametric") {
      ## --- parametric resampling: W^(b)_i ~ Mult(L_i, q_hat_i) -----------
      Wb <- matrix(0L, N, ncol(W))
      for (i in seq_len(N)) {
        if (L_obs[i] > 0L)
          Wb[i, ] <- as.integer(stats::rmultinom(1L, size = L_obs[i],
                                                 prob = q_hat[i, ]))
      }
      colnames(Wb) <- colnames(W)
      Xb <- X
    } else {
      ## --- non-parametric row resampling (legacy) ------------------------
      idx <- sample.int(N, N, replace = TRUE)
      Wb <- W[idx, , drop = FALSE]
      Xb <- X[idx, , drop = FALSE]
    }

    ctrl_b <- if (method == "parametric") control_boot else control
    ## (parametric: cold init for full first-stage uncertainty; see
    ## comment above the control_boot construction)
    fitb <- tryCatch(
      fit_gscamm(Wb, Xb, K = K, link = link,
                 gsca_space = fit$gsca_space %||% "alr",
                 gsca_ref = fit$gsca_ref %||% K,
                 init_phi = fit$init_phi %||% "kmeans",
                 polish = "none",   ## bootstrap only needs B; skip MAP polish
                 control = ctrl_b,
                 verbose = FALSE, seed = NULL),
      error = function(e) NULL
    )
    if (is.null(fitb)) return(NULL)
    ## align bootstrap fit to the original Phi
    perm <- align_components(fitb$Phi, fit$Phi)
    fitb_aligned <- .permute_fit(fitb, perm)

    eff_b <- tryCatch({
      if (method == "noise_augmented") {
        ## perturb the response (theta or R) in the ALR space, then map
        ## back to the simplex via softmax-with-ref before ALR-WLS.
        base_b <- if (use == "responsibilities" && !is.null(fitb_aligned$R))
          fitb_aligned$R else fitb_aligned$Theta
        base_b <- pmax(base_b, .Machine$double.eps)
        Nb <- nrow(base_b); Kb <- ncol(base_b)
        nr_b <- setdiff(seq_len(Kb), ref)
        alr_b <- log(base_b[, nr_b, drop = FALSE] / base_b[, ref])
        for (j in seq_along(nr_b)) {
          alr_b[, j] <- alr_b[, j] +
            stats::rnorm(Nb, sd = noise_sigma[j])
        }
        full <- matrix(0, Nb, Kb)
        full[, nr_b] <- alr_b
        m <- apply(full, 1, max); ez <- exp(full - m)
        theta_perturbed <- ez / rowSums(ez)
        fitb_perturbed <- fitb_aligned
        fitb_perturbed$Theta <- theta_perturbed
        fitb_perturbed$R     <- theta_perturbed
        covariate_effects(fitb_perturbed, ref = ref, level = level,
                          adjust = "none", use = "theta")
      } else {
        covariate_effects(fitb_aligned, ref = ref, level = level,
                          adjust = "none", use = use)
      }
    }, error = function(e) NULL)
    if (is.null(eff_b)) return(NULL)
    eff_b$B_alr
  }

  ## run replicates
  if (!is.null(seed)) set.seed(seed)
  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() %/% 2L)
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, library(gscamm))
      parallel::clusterExport(cl,
        varlist = c("W", "X", "N", "K", "P", "link", "control",
                    "control_boot", "fit", "ref", "level", "one_boot",
                    "method", "use", "noise_sigma",
                    "Theta_hat", "Phi_hat", "B_hat", "q_hat", "L_obs"),
        envir = environment())
      ## seed each worker
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)
      results <- parallel::parLapply(cl, seq_len(B), one_boot)
    } else {
      results <- parallel::mclapply(seq_len(B), one_boot, mc.cores = n_cores)
    }
  } else {
    results <- vector("list", B)
    for (b in seq_len(B)) {
      results[[b]] <- one_boot(b)
      if (verbose && b %% 10 == 0) message("bootstrap ", b, "/", B)
      if (!is.null(progress_callback)) progress_callback(b, B)
    }
  }

  ## stack into array, drop failures
  ok <- !vapply(results, is.null, logical(1))
  n_fail <- sum(!ok)
  if (n_fail > 0) message(n_fail, " bootstrap replicate(s) failed and were discarded.")
  results <- results[ok]
  if (length(results) < 10)
    stop("Too few successful bootstrap replicates (", length(results), ").")
  Bdraws <- array(NA_real_,
                  dim = c(length(results), P + 1L, K - 1L),
                  dimnames = list(NULL, cov_names, comp_names))
  for (b in seq_along(results)) Bdraws[b, , ] <- results[[b]]

  ## summarize
  alpha <- 1 - level
  qlo <- alpha / 2; qhi <- 1 - qlo

  est <- eff_orig$B_alr  ## point estimate from original fit
  boot_mean <- apply(Bdraws, c(2, 3), mean, na.rm = TRUE)
  bias_hat  <- boot_mean - est                        ## bootstrap bias estimate
  se        <- apply(Bdraws, c(2, 3), stats::sd, na.rm = TRUE)
  q_lo_emp  <- apply(Bdraws, c(2, 3), stats::quantile,
                     probs = qlo, na.rm = TRUE)
  q_hi_emp  <- apply(Bdraws, c(2, 3), stats::quantile,
                     probs = qhi, na.rm = TRUE)

  ## ---- compute CI matching the requested type ----------------------------
  ci <- switch(type,
    percentile = list(lo = q_lo_emp, hi = q_hi_emp),
    basic      = list(lo = 2 * est - q_hi_emp,
                      hi = 2 * est - q_lo_emp),
    bca = {
      ## bias-correction z0 from proportion of bootstrap below observed
      z0 <- array(NA_real_, dim = dim(est), dimnames = dimnames(est))
      for (j in seq_len(ncol(est))) for (i in seq_len(nrow(est))) {
        d <- Bdraws[, i, j]; d <- d[is.finite(d)]
        if (!length(d)) next
        p_lt <- mean(d < est[i, j])
        z0[i, j] <- stats::qnorm(min(max(p_lt, 1e-6), 1 - 1e-6))
      }
      ## acceleration via jackknife on the SECOND-STAGE WLS only (cheap)
      ## a = sum((Theta_bar_jk - Theta_jk)^3) / (6 * (sum((Theta_bar_jk - Theta_jk)^2))^(3/2))
      ## We approximate by jackknifing rows of (Theta, X) and recomputing
      ## ALR-WLS on the original fit, leaving the EM step untouched.
      a <- .bca_acceleration(fit, ref = ref, use = use)
      za <- stats::qnorm(c(qlo, qhi))
      lo <- est; hi <- est
      for (j in seq_len(ncol(est))) for (i in seq_len(nrow(est))) {
        if (!is.finite(z0[i, j])) {
          lo[i, j] <- q_lo_emp[i, j]; hi[i, j] <- q_hi_emp[i, j]; next
        }
        a_ij <- if (is.matrix(a)) a[i, j] else a
        adj <- function(zq) {
          v <- z0[i, j] + (z0[i, j] + zq) /
            (1 - a_ij * (z0[i, j] + zq))
          stats::pnorm(v)
        }
        plo <- adj(za[1]); phi <- adj(za[2])
        d <- Bdraws[, i, j]; d <- d[is.finite(d)]
        lo[i, j] <- stats::quantile(d, probs = plo, names = FALSE)
        hi[i, j] <- stats::quantile(d, probs = phi, names = FALSE)
      }
      list(lo = lo, hi = hi)
    }
  )
  ci_lo <- ci$lo; ci_hi <- ci$hi

  ## bias-corrected point estimate (and shifted CI), if requested
  if (bias_correct) {
    est_corr <- est - bias_hat
    shift <- est_corr - est
    ci_lo <- ci_lo + shift
    ci_hi <- ci_hi + shift
    est <- est_corr
  }

  ## percentile two-sided p-value: 2 * min(F*(0), 1 - F*(0))
  ## where F* is the bootstrap ECDF for each coefficient
  ## (kept on the un-shifted draws when bias_correct = TRUE; users who
  ## want adjusted p-values should rely on the CI excluding 0)
  pmat <- apply(Bdraws, c(2, 3), function(d) {
    d <- d[is.finite(d)]
    if (!length(d)) return(NA_real_)
    f0 <- mean(d <= 0)
    2 * min(f0, 1 - f0)
  })
  ## clamp p-values away from exact 0 / 1 to a Monte-Carlo lower bound
  Bn <- length(results)
  pmin_mc <- 1 / (Bn + 1)
  pmat[is.finite(pmat) & pmat < pmin_mc] <- pmin_mc

  ## flatten to long data.frame consistent with covariate_effects()
  K1 <- K - 1L
  long <- data.frame(
    component  = rep(comp_names, each = P + 1L),
    covariate  = rep(cov_names, times = K1),
    estimate   = as.numeric(est),
    std.error  = as.numeric(se),
    statistic  = as.numeric(est) / pmax(as.numeric(se), .Machine$double.eps),
    p.value    = as.numeric(pmat),
    conf.low   = as.numeric(ci_lo),
    conf.high  = as.numeric(ci_hi),
    odds.ratio = exp(as.numeric(est)),
    or.low     = exp(as.numeric(ci_lo)),
    or.high    = exp(as.numeric(ci_hi)),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
  is_int <- long$covariate == "(Intercept)"
  long$p.adj <- NA_real_
  long$p.adj[!is_int] <- stats::p.adjust(long$p.value[!is_int], method = adjust)
  long <- long[, c("component", "covariate", "estimate", "std.error",
                   "statistic", "p.value", "p.adj", "conf.low",
                   "conf.high", "odds.ratio", "or.low", "or.high")]

  out <- list(
    coefficients = long,
    B_alr = est,
    vcov = NULL,
    B_draws = Bdraws,
    B_draws_failed = n_fail,
    bias = bias_hat,
    method = paste0("bootstrap_", method, "_", type,
                    if (bias_correct) "_bc" else ""),
    resample_method = method,
    type = type, bias_correct = bias_correct,
    ref = ref, level = level, adjust = adjust,
    K = K, P = P, B = length(results)
  )
  class(out) <- c("gscamm_effects_boot", "gscamm_effects")
  out
}

## -- BCa acceleration via leave-one-out jackknife on the second-stage
## ALR-WLS estimator, using the original fit's Theta. Cheap because it
## skips re-running the EM-GSCA loop for each jackknife replicate.
.bca_acceleration <- function(fit, ref, use = "theta") {
  N <- fit$N; K <- fit$K; P <- fit$P
  Theta_full <- if (use == "responsibilities" && !is.null(fit$R))
    fit$R else fit$Theta
  Theta_full <- pmax(Theta_full, .Machine$double.eps)
  X <- if (!is.null(fit$X_std)) fit$X_std else fit$X
  D <- cbind(`(Intercept)` = 1, X)
  non_ref <- setdiff(seq_len(K), ref)
  B_jk <- array(NA_real_, c(N, P + 1L, K - 1L))

  for (i in seq_len(N)) {
    Di <- D[-i, , drop = FALSE]
    Th <- Theta_full[-i, , drop = FALSE]
    for (j in seq_along(non_ref)) {
      k <- non_ref[j]
      y <- log(Th[, k] / Th[, ref])
      w <- Th[, k] * Th[, ref]
      f <- stats::lm.wfit(Di, y, w)
      B_jk[i, , j] <- f$coefficients
    }
  }
  Bbar <- apply(B_jk, c(2, 3), mean, na.rm = TRUE)
  num <- apply(sweep(B_jk, c(2, 3), Bbar, "-")^3, c(2, 3), sum, na.rm = TRUE)
  den <- 6 * apply(sweep(B_jk, c(2, 3), Bbar, "-")^2,
                   c(2, 3), sum, na.rm = TRUE)^(3 / 2)
  out <- ifelse(den > 0, num / den, 0)
  out
}
