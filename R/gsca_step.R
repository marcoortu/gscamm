## ---------------------------------------------------------------------------
## GSCA step: update path coefficients B and component scores Gamma.
##
## Reference: Ortu and Frigau (2026), Section 2.3, equations (8)-(9), and
## Remark 1 (geometric ridge schedule).
##
## Identifiability convention (paper Section 2.2): the path coefficient
## matrix is parametrized as B = [B_{-K}, 0], i.e. the reference component
## carries a structural zero column. Only the non-reference block
## B_{-K} \in R^{P x (K-1)} is estimated; the GSCA step never updates the
## reference column. This is enforced by both the simplex- and ALR-space
## update routines below.
## ---------------------------------------------------------------------------

#' Ridge update for the simplex-space path coefficients
#'
#' Regresses the non-reference responsibilities r_{ik} on the standardized
#' covariates with a ridge penalty. The reference column is structurally
#' zero (identifiability constraint, paper Section 2.2).
#'
#' @param X_std standardized covariate matrix \eqn{N \times P}.
#' @param R responsibilities matrix \eqn{N \times K}.
#' @param lambda_B current ridge parameter (>= 0).
#' @param ref reference component (default last).
#' @param XtX precomputed \eqn{X^\top X} (optional cache).
#'
#' @return numeric matrix \eqn{P \times K} with a zero column at
#'   position \code{ref}.
#' @keywords internal
.gsca_update_B <- function(X_std, R, lambda_B, ref = ncol(R), XtX = NULL) {
  P <- ncol(X_std); K <- ncol(R)
  if (is.null(XtX)) XtX <- crossprod(X_std)
  A <- XtX + lambda_B * diag(P)
  rhs <- crossprod(X_std, R[, -ref, drop = FALSE])     ## P x (K-1)
  B_minus <- solve(A, rhs)
  B <- matrix(0, P, K)
  B[, -ref] <- B_minus
  B
}

#' Ridge update for the path coefficients in ALR space
#'
#' GSCA step formulated as a penalized log-ratio projection of the
#' posterior responsibilities onto the covariate space (paper Section 2.3,
#' Remark on the ALR-coordinate update). Each non-reference responsibility
#' is regularized by a small constant \eqn{\delta} before the additive
#' log-ratio transform:
#' \deqn{\tilde r_{ik} = (r_{ik} + \delta) / \sum_\ell (r_{i\ell} + \delta).}
#' The non-reference block of \eqn{B} is then the ridge multivariate
#' regression of the ALR responses on the standardized covariates.
#'
#' Returns a \eqn{P \times K} matrix with a zero column at the reference
#' position (identifiability constraint).
#'
#' @param X_std standardized covariate matrix \eqn{N \times P}.
#' @param R responsibilities matrix \eqn{N \times K}.
#' @param lambda_B current ridge parameter (>= 0).
#' @param ref reference component (default last).
#' @param delta non-negative ALR regularization constant; must be positive
#'   when responsibilities can take exact zeros (default \code{1e-2}).
#' @param XtX precomputed \eqn{X^\top X} (optional).
#' @keywords internal
.gsca_update_B_alr <- function(X_std, R, lambda_B, ref = ncol(R),
                               delta = 1e-2, XtX = NULL) {
  P <- ncol(X_std); K <- ncol(R)
  if (is.null(XtX)) XtX <- crossprod(X_std)
  ## delta-regularize the responsibilities, then ALR-transform with ref
  ## as denominator (paper roadmap point 3)
  Rs <- R + delta
  Rs <- Rs / rowSums(Rs)
  ALR <- log(Rs[, -ref, drop = FALSE] / Rs[, ref])
  A <- XtX + lambda_B * diag(P)
  B_alr <- solve(A, crossprod(X_std, ALR))             ## P x (K-1)
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
