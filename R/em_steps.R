## ---------------------------------------------------------------------------
## EM updates for the mixture-of-multinomials part of GSCA-MM.
##
## Reference: Ortu and Frigau (2026), Section 2.3 "Penalized Likelihood
## Estimation". The E-step computes document-level posterior responsibilities
## (Equation 6, aggregated form). The M-step updates the component-category
## distributions Phi with Dirichlet smoothing alpha (Equation 7).
##
## All operations are done in log-space to avoid underflow when V is large
## or document lengths are heterogeneous.
## ---------------------------------------------------------------------------

#' E-step: aggregated document-level responsibilities
#'
#' Computes \eqn{r_{ik} \propto \theta_{ik} \prod_v \phi_{kv}^{n_{iv}}},
#' normalized rowwise to sum to one. Uses the log-sum-exp trick.
#'
#' @param W non-negative numeric matrix \eqn{N \times V} of counts.
#' @param Theta numeric matrix \eqn{N \times K}.
#' @param Phi numeric matrix \eqn{K \times V}.
#'
#' @return numeric matrix \eqn{N \times K} of responsibilities.
#' @keywords internal
.e_step <- function(W, Theta, Phi) {
  ## log Phi with floor for numerical stability
  log_phi <- log(pmax(Phi, .Machine$double.xmin))
  log_theta <- log(pmax(Theta, .Machine$double.xmin))

  ## sum_v n_iv * log phi_kv  ==  W %*% t(log_phi)  -> N x K
  log_rk <- log_theta + W %*% t(log_phi)

  ## normalize each row by log-sum-exp
  m <- apply(log_rk, 1, max)
  z <- exp(log_rk - m)
  z / rowSums(z)
}

#' M-step: update of component-category distributions Phi
#'
#' \eqn{\phi_{kv}^{(t+1)} \propto \sum_i r_{ik} n_{iv} + \alpha}
#' (Equation 7), with row-normalization.
#'
#' @param W non-negative numeric matrix \eqn{N \times V} of counts.
#' @param R responsibilities matrix \eqn{N \times K}.
#' @param alpha non-negative Dirichlet smoothing constant.
#'
#' @return numeric matrix \eqn{K \times V} on the simplex.
#' @keywords internal
.m_step_phi <- function(W, R, alpha = 0.01) {
  ## sufficient statistics: t(R) %*% W -> K x V
  num <- crossprod(R, W) + alpha
  num / rowSums(num)
}
