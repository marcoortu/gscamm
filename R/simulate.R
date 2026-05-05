## ---------------------------------------------------------------------------
## Synthetic data generator implementing the simulation design of Section 4
## of Ortu and Frigau (2026), generalized to mixture-of-multinomials models.
##
## Note: the GSCA-MM *model* used by `fit_gscamm` is deterministic
## (eta = X B_{-K}, theta = inverse-ALR(eta); see fit_gscamm.R). The
## additive Gaussian term in the ALR predictor below is part of the
## *data-generating mechanism* only -- it perturbs the simulated theta
## away from a degenerate covariate-driven manifold so that we can study
## recovery of B from noisy mixture weights, and is the analogue of the
## first-stage uncertainty addressed by the bootstrap and by the
## responsibility-based plug-in inference at fit time.
##
## Three scenarios are supported:
##  - "baseline":       logistic-normal, sigma = 0.3, lengths Poisson(1000)
##  - "high_covariate": Dirichlet(alpha * theta) with alpha = 50
##  - "high_sparsity":  ALR-then-threshold at 0.03, lengths Poisson(20)
## ---------------------------------------------------------------------------

#' Simulate a GSCA-MM dataset
#'
#' Generates a count matrix \eqn{W} together with a covariate matrix
#' \eqn{X} and the underlying ground-truth parameters, mirroring the
#' simulation design of Section 4 of the paper.
#'
#' Covariate effects are encoded through an ALR parameterization: for
#' each observation a vector of \eqn{K-1} log-ratios is drawn as
#' \eqn{\eta_i = X_i^\top \beta + \varepsilon_i}, augmented with a zero
#' reference component, and mapped to mixture proportions through the
#' inverse ALR transform (softmax with the last component fixed to zero).
#'
#' @param N number of observations.
#' @param V number of categories.
#' @param K number of components.
#' @param P number of covariates.
#' @param scenario one of \code{"baseline"} (default),
#'   \code{"high_covariate"}, \code{"high_sparsity"}.
#' @param sigma residual sd of the ALR linear predictor (default 0.3).
#' @param dirichlet_alpha concentration parameter for the
#'   \code{"high_covariate"} scenario (default 50).
#' @param sparsity_thr threshold under which proportions are zeroed in the
#'   \code{"high_sparsity"} scenario (default 0.03).
#' @param doc_length_mean Poisson mean for observation lengths.
#' @param doc_length_min lower truncation for observation lengths
#'   (default 5).
#' @param phi_alpha Dirichlet concentration for sampling
#'   component-category distributions (default 0.1).
#' @param beta_scale standard deviation of the Gaussian draw used to
#'   generate the true \eqn{\beta} matrix (default 1).
#' @param seed optional integer seed.
#'
#' @return a list with elements:
#'   \item{W}{N x V count matrix.}
#'   \item{X}{N x P covariate matrix.}
#'   \item{Theta}{N x K true mixture proportions.}
#'   \item{Phi}{K x V true component-category matrix.}
#'   \item{beta}{P x (K-1) true ALR coefficient matrix (reference = K).}
#'   \item{lengths}{N integer vector of observation lengths.}
#'   \item{scenario}{the scenario name.}
#'
#' @export
simulate_gscamm <- function(N = 1000, V = 500, K = 10, P = 8,
                            scenario = c("baseline", "high_covariate",
                                         "high_sparsity"),
                            sigma = 0.3,
                            dirichlet_alpha = 50,
                            sparsity_thr = 0.03,
                            doc_length_mean = NULL,
                            doc_length_min = 5,
                            phi_alpha = 0.1,
                            beta_scale = 1,
                            seed = NULL) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(doc_length_mean))
    doc_length_mean <- if (scenario == "high_sparsity") 20 else 1000

  ## --- covariates ----------------------------------------------------------
  X <- matrix(stats::rnorm(N * P), N, P)
  X_std <- .standardize(X)

  ## --- true ALR coefficients beta (P x (K-1)) ------------------------------
  beta <- matrix(stats::rnorm(P * (K - 1L), sd = beta_scale), P, K - 1L)

  ## --- ALR linear predictor and back-transform to Theta --------------------
  eta <- X_std %*% beta + matrix(stats::rnorm(N * (K - 1L), sd = sigma),
                                 N, K - 1L)
  ## inverse ALR with last component as reference
  eta_full <- cbind(eta, 0)
  m <- apply(eta_full, 1, max)
  Theta <- exp(eta_full - m)
  Theta <- Theta / rowSums(Theta)

  if (scenario == "high_covariate") {
    Theta <- t(apply(Theta, 1, function(p) {
      a <- dirichlet_alpha * pmax(p, 1e-12)
      g <- stats::rgamma(length(a), a, 1)
      g / sum(g)
    }))
  } else if (scenario == "high_sparsity") {
    Theta[Theta < sparsity_thr] <- 0
    rs <- rowSums(Theta)
    rs[rs == 0] <- 1
    Theta <- Theta / rs
  }

  ## --- component-category distributions Phi --------------------------------
  Phi <- matrix(stats::rgamma(K * V, phi_alpha, 1), K, V)
  Phi <- Phi / rowSums(Phi)

  ## --- sample observations from the multinomial mixture --------------------
  q <- Theta %*% Phi
  L <- pmax(stats::rpois(N, doc_length_mean), doc_length_min)
  W <- matrix(0L, N, V)
  for (i in seq_len(N)) {
    W[i, ] <- as.integer(stats::rmultinom(1, L[i], prob = q[i, ]))
  }

  list(W = W, X = X, Theta = Theta, Phi = Phi, beta = beta,
       lengths = L, scenario = scenario)
}
