## ---------------------------------------------------------------------------
## Post-estimation: ALR-WLS regression of mixture proportions on covariates.
##
## Reference: Ortu and Frigau (2026), Section 2.4 "Estimating Covariate
## Effects on Mixture Proportions". For each non-reference component k,
## fit a weighted linear model
##
##     y_ik = beta_0k + x_i' beta_k + eps_ik
##
## with y_ik = log(theta_ik / theta_iK) and weights w_ik = theta_ik *
## theta_iK. Sandwich/SE-from-WLS, Wald p-values, BH adjustment, and
## odds-ratio CIs follow standard practice. The two-stage inference scheme
## with consistency / asymptotic normality is established in
## Theorem 3 (Section 3.3) of the paper.
## ---------------------------------------------------------------------------

#' Estimate covariate effects on mixture proportions
#'
#' Performs the additive-log-ratio (ALR) weighted least-squares regression
#' of fitted mixture proportions on covariates, as described in Section
#' 2.4 of the GSCA-MM paper. For each non-reference component
#' \eqn{k = 1, \ldots, K-1} a weighted linear model is fitted with
#' response \eqn{y_{ik} = \log(\hat{\theta}_{ik}/\hat{\theta}_{iK})} and
#' weights \eqn{w_{ik} = \hat{\theta}_{ik} \hat{\theta}_{iK}}, where
#' component \code{ref} plays the role of the reference category.
#'
#' @param object a \code{gscamm} object returned by \code{\link{fit_gscamm}}.
#' @param ref reference component index (default the last one, \eqn{K}).
#' @param level confidence level for intervals on coefficients and odds
#'   ratios (default 0.95).
#' @param adjust multiple-testing adjustment method passed to
#'   \code{\link[stats]{p.adjust}} (default \code{"BH"}).
#' @param standardize logical: if \code{TRUE} (default) the regression uses
#'   the standardized covariate matrix that the model was fit on, so
#'   coefficients are on the standardized scale. Setting \code{FALSE}
#'   re-fits the regression on the original (un-standardized) covariates
#'   if they are available in \code{object$X_raw}.
#' @param X_new optional alternative design matrix to use instead of the
#'   one stored in the fit. Must have the same number of rows as the data.
#' @param use response basis for the ALR-WLS *point* estimate. Either
#'   \code{"theta"} (default; the converged mixture proportions, matching
#'   Section 2.4 of the paper) or \code{"responsibilities"} (the
#'   unit-level posterior responsibilities from the final E-step). Theta
#'   is concentrated by the deterministic link and is the canonical
#'   target of the structural ALR coefficients, so it gives the point
#'   estimate with smallest bias against the data-generating beta;
#'   responsibilities give a noisier estimate of similar mean but with
#'   different variance properties.
#' @param var_method estimator for the second-stage variance: either
#'   \code{"hybrid"} (default) or \code{"hc0"}. The hybrid estimator
#'   computes \eqn{\hat\sigma_k^2} from R-based residuals
#'   \eqn{e_i = \log(r_{ik}/r_{i,\mathrm{ref}}) - x_i^\top \hat\beta_k}
#'   and forms \eqn{V = \hat\sigma_k^2 (D^\top W D)^{-1}}. This is the
#'   recommended workflow under the deterministic link: Theta yields a
#'   small-bias point estimate while R supplies the first-stage residual
#'   variance that is otherwise lost. The HC0 sandwich is retained for
#'   parity with the original paper.
#'
#' @return an object of class \code{gscamm_effects}, a list with components:
#'   \item{coefficients}{long-format data frame with columns
#'     \code{component}, \code{covariate}, \code{estimate}, \code{std.error},
#'     \code{statistic}, \code{p.value}, \code{p.adj}, \code{conf.low},
#'     \code{conf.high}, \code{odds.ratio}, \code{or.low}, \code{or.high}.}
#'   \item{B_alr}{\eqn{(P+1) \times (K-1)} matrix of estimated coefficients
#'     (intercept in row 1).}
#'   \item{vcov}{list of (P+1)x(P+1) sandwich covariance matrices, one per
#'     non-reference component.}
#'   \item{ref}{the reference component index.}
#'   \item{level, adjust}{settings used.}
#'
#' @references Ortu and Frigau (2026), Section 2.4 and Theorem 3.
#' @export
covariate_effects <- function(object,
                              ref = object$K,
                              level = 0.95,
                              adjust = c("BH", "holm", "hochberg",
                                         "hommel", "bonferroni", "BY",
                                         "fdr", "none"),
                              standardize = TRUE,
                              X_new = NULL,
                              use = c("theta", "responsibilities"),
                              var_method = c("hybrid", "hc0")) {
  if (!inherits(object, "gscamm"))
    stop("`object` must be a gscamm fit.")
  adjust     <- match.arg(adjust)
  use        <- match.arg(use)
  var_method <- match.arg(var_method)
  K <- object$K
  if (ref < 1 || ref > K) stop("`ref` out of range.")

  ## response basis for the ALR-WLS point estimate
  base <- switch(use,
                 theta            = object$Theta,
                 responsibilities = object$R %||% object$Theta)
  Theta <- pmax(base, .Machine$double.eps)
  ## responsibility basis for the hybrid first-stage residual variance
  R_base <- pmax(object$R %||% object$Theta, .Machine$double.eps)

  ## design matrix for the second-stage regression
  if (!is.null(X_new)) {
    Xd <- as.matrix(X_new)
    if (nrow(Xd) != object$N)
      stop("X_new must have the same number of rows as the fitted data.")
  } else if (standardize) {
    Xd <- object$X_std
  } else if (!is.null(object$X_raw)) {
    Xd <- object$X_raw
  } else {
    Xd <- object$X_std
  }
  ## column names default to V1..VP
  if (is.null(colnames(Xd)))
    colnames(Xd) <- paste0("X", seq_len(ncol(Xd)))

  P <- ncol(Xd)
  D <- cbind(`(Intercept)` = 1, Xd)  ## N x (P+1)
  cov_names <- colnames(D)

  non_ref <- setdiff(seq_len(K), ref)
  K1 <- length(non_ref)

  B_alr <- matrix(NA_real_, P + 1L, K1,
                  dimnames = list(cov_names, paste0("comp", non_ref)))
  vcov_list <- vector("list", K1)
  names(vcov_list) <- paste0("comp", non_ref)

  long_rows <- vector("list", K1)
  z_crit <- stats::qnorm(1 - (1 - level) / 2)

  for (j in seq_along(non_ref)) {
    k <- non_ref[j]
    y <- log(Theta[, k] / Theta[, ref])
    w <- Theta[, k] * Theta[, ref]
    ## weighted least squares via lm.wfit
    fit_k <- stats::lm.wfit(D, y, w)
    bhat <- fit_k$coefficients
    bhat[is.na(bhat)] <- 0
    B_alr[, j] <- bhat

    XtWX <- crossprod(D, w * D)
    XtWX_inv <- tryCatch(solve(XtWX), error = function(...) MASS_ginv(XtWX))

    if (var_method == "hc0") {
      ## sandwich-style WLS variance:
      ##   Var(beta) = (X' W X)^(-1) X' W diag(e^2) W X (X' W X)^(-1)
      e <- as.numeric(fit_k$residuals)
      meat <- crossprod(D, (w^2 * e^2) * D)
      Vmat <- XtWX_inv %*% meat %*% XtWX_inv
    } else {
      ## hybrid estimator: point estimate from Theta-based WLS, residual
      ## variance from responsibility-based residuals (paper Remark on
      ## first-stage uncertainty). Captures the variability that the
      ## deterministic link otherwise discards.
      y_R <- log(R_base[, k] / R_base[, ref])
      e_R <- as.numeric(y_R - D %*% bhat)
      df_k <- max(length(e_R) - ncol(D), 1L)
      sigma2_k <- sum(w * e_R^2) / df_k
      Vmat <- sigma2_k * XtWX_inv
    }
    dimnames(Vmat) <- list(cov_names, cov_names)
    vcov_list[[j]] <- Vmat

    se <- sqrt(pmax(diag(Vmat), 0))
    stat <- bhat / se
    p <- 2 * stats::pnorm(-abs(stat))
    ci_lo <- bhat - z_crit * se
    ci_hi <- bhat + z_crit * se

    long_rows[[j]] <- data.frame(
      component  = paste0("comp", k),
      covariate  = cov_names,
      estimate   = bhat,
      std.error  = se,
      statistic  = stat,
      p.value    = p,
      conf.low   = ci_lo,
      conf.high  = ci_hi,
      odds.ratio = exp(bhat),
      or.low     = exp(ci_lo),
      or.high    = exp(ci_hi),
      row.names  = NULL,
      stringsAsFactors = FALSE
    )
  }

  coefs <- do.call(rbind, long_rows)
  ## adjust p-values across all (covariate, component) pairs jointly,
  ## excluding the intercept (consistent with Equation (12) of the paper)
  is_int <- coefs$covariate == "(Intercept)"
  p_adj <- rep(NA_real_, nrow(coefs))
  p_adj[!is_int] <- stats::p.adjust(coefs$p.value[!is_int], method = adjust)
  coefs$p.adj <- p_adj
  coefs <- coefs[, c("component", "covariate", "estimate", "std.error",
                     "statistic", "p.value", "p.adj", "conf.low",
                     "conf.high", "odds.ratio", "or.low", "or.high")]

  out <- list(
    coefficients = coefs,
    B_alr = B_alr,
    vcov = vcov_list,
    ref = ref,
    level = level,
    adjust = adjust,
    K = K, P = P
  )
  class(out) <- "gscamm_effects"
  out
}

## tiny pseudo-inverse fallback to avoid hard MASS dependency
MASS_ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  d <- s$d
  d_inv <- ifelse(d > tol * max(d), 1 / d, 0)
  s$v %*% (d_inv * t(s$u))
}
