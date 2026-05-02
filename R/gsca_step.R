## ---------------------------------------------------------------------------
## GSCA step: update path coefficients B and component scores Gamma.
##
## Reference: Ortu and Frigau (2026), Section 2.3, equations (8)-(9), and
## Remark 1 (geometric ridge schedule). For each component k, B_k is the
## ridge-penalized least-squares fit of the responsibilities r_k on the
## standardized covariates X_std:
##
##     hat b_k = (X' X + lambda_B I)^(-1) X' r_k.
##
## Gamma is then the row-binding of X_std and the component scores
## X_std %*% B.
## ---------------------------------------------------------------------------

#' Ridge update for the path coefficient matrix B
#'
#' @param X_std standardized covariate matrix \eqn{N \times P}.
#' @param R responsibilities matrix \eqn{N \times K}.
#' @param lambda_B current ridge parameter (>= 0).
#' @param XtX precomputed \eqn{X^\top X} (optional cache); recomputed if
#'   \code{NULL}.
#'
#' @return numeric matrix \eqn{P \times K} of path coefficients.
#' @keywords internal
.gsca_update_B <- function(X_std, R, lambda_B, XtX = NULL) {
  P <- ncol(X_std)
  if (is.null(XtX)) XtX <- crossprod(X_std)
  A <- XtX + lambda_B * diag(P)
  ## solve once for the right-hand side X' R (P x K)
  rhs <- crossprod(X_std, R)
  solve(A, rhs)
}

#' Ridge update for the path coefficient matrix B in ALR space
#'
#' Alternative GSCA step that regresses the additive log-ratios of the
#' responsibilities (with a reference component) on the standardized
#' covariates. Compared with the simplex-space update of equation (8),
#' this variant aligns the linear regression with the natural geometry of
#' compositional data and removes the structural bias induced by
#' regressing simplex-valued responsibilities on covariates.
#'
#' Returns a \eqn{P \times K} matrix obtained by re-injecting a column of
#' zeros for the reference component, so that downstream code sees the
#' same shape as in the simplex-space variant.
#'
#' @param X_std standardized covariate matrix \eqn{N \times P}.
#' @param R responsibilities matrix \eqn{N \times K}.
#' @param lambda_B current ridge parameter (>= 0).
#' @param ref reference component (default last).
#' @param XtX precomputed \eqn{X^\top X} (optional).
#' @keywords internal
.gsca_update_B_alr <- function(X_std, R, lambda_B, ref = ncol(R),
                               XtX = NULL) {
  P <- ncol(X_std); K <- ncol(R)
  if (is.null(XtX)) XtX <- crossprod(X_std)
  Rs <- pmax(R, .Machine$double.eps)
  ## ALR-transform with ref as denominator
  ALR <- log(Rs[, -ref, drop = FALSE] / Rs[, ref])
  A <- XtX + lambda_B * diag(P)
  B_alr <- solve(A, crossprod(X_std, ALR))             ## P x (K-1)
  ## inject zeros for the reference column to preserve shape P x K
  B <- matrix(0, P, K)
  B[, -ref] <- B_alr
  B
}

#' Build the component-score matrix Gamma
#'
#' Equation (9): \eqn{\boldsymbol{\Gamma} = [X_{std}, X_{std} B]}, an
#' \eqn{N \times (P+K)} matrix whose first \eqn{P} columns are the
#' standardized covariates and last \eqn{K} columns are the component
#' scores.
#'
#' @param X_std standardized covariate matrix \eqn{N \times P}.
#' @param B path coefficient matrix \eqn{P \times K}.
#'
#' @return numeric matrix \eqn{N \times (P+K)}.
#' @keywords internal
.build_gamma <- function(X_std, B) {
  cbind(X_std, X_std %*% B)
}

#' Geometric ridge schedule \eqn{\lambda_B^{(t)} = \lambda_0 \rho^t}
#'
#' Remark 1 of the paper. Defaults are \eqn{\lambda_0 = 1}, \eqn{\rho = 0.9}.
#'
#' @param t iteration index (zero-based).
#' @param lambda0 initial ridge parameter.
#' @param rho decay rate in (0, 1).
#'
#' @return positive scalar.
#' @keywords internal
.lambda_schedule <- function(t, lambda0 = 1, rho = 0.9) {
  lambda0 * rho^t
}
