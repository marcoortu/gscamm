## ---------------------------------------------------------------------------
## Per-observation MAP polish in ALR space.
##
## After EM-GSCA convergence the structural estimator
##   Theta_i = inverse-ALR(X_i' B_{-K})
## is X-only and deterministic, so it does not capture the row-level
## stochastic component implicit in the data-generating mechanism (paper
## Section 4.1, eq. 13: eta_i = X_i' beta + epsilon_i). The token-level
## responsibilities R = E[Z|...] over-peak with informative documents and
## are also a poor estimator of theta_i.
##
## The MAP estimate combines both signals: it solves, per row i,
##   eta_i^MAP = argmax  L(eta; W_i, Phi) - 0.5 (eta - mu_i)' Sigma^-1 (eta - mu_i)
## with mu_i = X_i' B_{-K} and Sigma = diag(sigma^2_k) estimated from ALR
## residuals of R against the structural mean. The resulting
##   Theta_map_i = inverse-ALR(eta_map_i)
## is the apples-to-apples counterpart of LDA gamma and STM theta:
## covariate-aware prior + data-informed likelihood.
##
## Used as the default theta estimator for gscamm in the simulation study
## (see fit_gscamm(polish = "map")), without affecting the structural B
## estimate or its sandwich variance.
## ---------------------------------------------------------------------------

#' MAP polish of mixture proportions in ALR space
#'
#' For each observation, solves a Newton optimization to combine the
#' covariate-driven structural mean
#' \eqn{\boldsymbol{\mu}_i = X_i^\top \widehat{B}_{-K}} with the
#' mixture-of-multinomials likelihood under \eqn{\widehat{\Phi}} and a
#' Gaussian ALR prior with diagonal variance \code{sigma2}. Returns the
#' MAP score matrix and its inverse-ALR mixture proportions.
#'
#' Optimization uses Newton steps with a Fisher-information approximation
#' to the Hessian and a backtracking line search. Cost is
#' \eqn{O(N\,K^2\,T)} with \eqn{T \le} \code{max_iter}; typical \eqn{T \le 10}.
#'
#' @param W non-negative count matrix \eqn{N \times V}.
#' @param Phi component-category matrix \eqn{K \times V} on the simplex.
#' @param mu prior mean in ALR space, \eqn{N \times (K-1)}.
#' @param sigma2 length-\eqn{K-1} positive prior variances.
#' @param ref reference component index (default \code{nrow(Phi)}).
#' @param max_iter Newton iterations per row (default 20).
#' @param tol max-absolute-change tolerance on \eqn{\eta} (default 1e-6).
#'
#' @return list with components:
#'   \item{eta}{\eqn{N \times (K-1)} MAP scores in ALR space.}
#'   \item{Theta}{\eqn{N \times K} mixture proportions = inverse-ALR(eta).}
#'   \item{converged}{logical vector of length \code{N}.}
#'   \item{iters}{integer vector of length \code{N}.}
#'
#' @keywords internal
.eta_map_polish <- function(W, Phi, mu, sigma2, ref = nrow(Phi),
                            max_iter = 20L, tol = 1e-6) {
  W <- as.matrix(W)
  Phi <- as.matrix(Phi)
  mu <- as.matrix(mu)
  N <- nrow(W); V <- ncol(W); K <- nrow(Phi)
  Km1 <- K - 1L
  if (ncol(mu) != Km1)
    stop(".eta_map_polish: mu must have K-1 columns.")
  if (length(sigma2) != Km1)
    stop(".eta_map_polish: sigma2 must have length K-1.")
  if (any(sigma2 <= 0))
    stop(".eta_map_polish: sigma2 entries must be positive.")
  if (ref < 1L || ref > K)
    stop(".eta_map_polish: ref out of range.")

  nr   <- setdiff(seq_len(K), ref)             ## non-reference indices in 1..K
  prec <- 1 / sigma2                           ## length Km1, prior precisions

  eta   <- mu                                  ## N x Km1, init at prior mean
  iters <- integer(N)
  conv  <- logical(N)

  L_i   <- rowSums(W)
  empty <- L_i == 0
  if (any(empty)) {
    conv[empty]  <- TRUE
    iters[empty] <- 0L
  }

  for (i in seq_len(N)[!empty]) {
    et   <- as.numeric(eta[i, ])
    mu_i <- as.numeric(mu[i, ])
    W_i  <- as.numeric(W[i, ])
    Li   <- L_i[i]

    J_cur <- .map_obj_row(et, mu_i, prec, W_i, Phi, nr, ref, K)

    for (t in seq_len(max_iter)) {
      ## --- forward pass: theta, q, gradient, Fisher info -------------------
      eta_full <- numeric(K); eta_full[nr] <- et   ## eta_ref = 0 by convention
      m <- max(eta_full)
      th <- exp(eta_full - m); th <- th / sum(th)  ## length K

      q  <- as.numeric(crossprod(th, Phi))         ## length V; q_v = th' Phi_,v
      q  <- pmax(q, .Machine$double.xmin)

      Wq <- W_i / q                                ## length V

      ## gradient of log-lik wrt eta_l (l in nr):
      ##   g_l = th_l * ( sum_v W_iv * Phi_lv / q_v  -  L_i )
      A <- as.numeric(Phi %*% Wq)                  ## length K
      g_full <- th * (A - Li)                      ## length K (g_ref is dropped)
      g <- g_full[nr]                              ## length Km1

      ## Fisher information block (Km1 x Km1):
      ##   I_{lm} = L_i * th_l * th_m * sum_v (Phi_lv - q_v)(Phi_mv - q_v) / q_v
      E   <- Phi - matrix(q, K, V, byrow = TRUE)   ## K x V
      Esc <- E / matrix(sqrt(q), K, V, byrow = TRUE)
      M   <- tcrossprod(Esc)                       ## K x K
      I_full <- Li * (th %o% th) * M               ## K x K
      I_red  <- I_full[nr, nr, drop = FALSE]

      ## --- Newton step on the MAP system (I_red + diag(prec)) d = g - prec*(eta - mu)
      H_pos <- I_red + diag(prec, Km1)
      rhs   <- g - prec * (et - mu_i)
      delta <- tryCatch(
        solve(H_pos, rhs),
        error = function(e) solve(H_pos + diag(1e-6, Km1), rhs)
      )

      ## --- backtracking: ensure J non-decreasing
      step  <- 1.0
      moved <- FALSE
      for (bt in 1:8) {
        et_new <- et + step * delta
        J_new  <- .map_obj_row(et_new, mu_i, prec, W_i, Phi, nr, ref, K)
        if (is.finite(J_new) && J_new >= J_cur - 1e-10) {
          et    <- et_new
          J_cur <- J_new
          moved <- TRUE
          break
        }
        step <- step / 2
      }

      iters[i] <- t
      if (!moved) break                            ## stuck: keep current et
      if (max(abs(step * delta)) < tol) {
        conv[i] <- TRUE
        break
      }
    }
    eta[i, ] <- et
  }

  ## --- inverse-ALR back to the simplex ---------------------------------------
  eta_full_mat        <- matrix(0, N, K)
  eta_full_mat[, nr]  <- eta
  m  <- apply(eta_full_mat, 1, max)
  Th <- exp(eta_full_mat - m)
  Th <- Th / rowSums(Th)

  list(eta = eta, Theta = Th, converged = conv, iters = iters)
}

#' MAP objective evaluated at a single row's eta (log-lik minus prior penalty)
#' @keywords internal
.map_obj_row <- function(et, mu_i, prec, W_i, Phi, nr, ref, K) {
  eta_full <- numeric(K); eta_full[nr] <- et
  m  <- max(eta_full)
  th <- exp(eta_full - m); th <- th / sum(th)
  q  <- as.numeric(crossprod(th, Phi))
  q  <- pmax(q, .Machine$double.xmin)
  ll  <- sum(W_i * log(q))
  pen <- 0.5 * sum(prec * (et - mu_i)^2)
  ll - pen
}

#' Estimate per-component prior variance from ALR residuals of responsibilities
#'
#' Given the converged responsibilities \eqn{R} and the structural mean
#' \eqn{X_i^\top \widehat{B}_{-K}}, computes
#' \deqn{\widehat\sigma^2_k = \frac{1}{N - P} \sum_i \big(\mathrm{alr}_k(\tilde R_i) - X_i^\top \widehat{B}_{-K,k}\big)^2,}
#' where \eqn{\tilde R_i} is \code{R} after \code{delta}-smoothing and ALR
#' is taken with reference \code{ref}. A small floor avoids degeneracies.
#'
#' @param R responsibilities matrix \eqn{N \times K}.
#' @param X_std standardized covariates \eqn{N \times P}.
#' @param B_minus structural coefficients \eqn{P \times (K-1)} (non-reference).
#' @param ref reference component index.
#' @param delta smoothing constant for ALR (matches control$delta).
#' @param P_eff effective number of estimated covariates (defaults to ncol(X_std)).
#' @param floor minimum allowed variance (default 1e-8).
#'
#' @return numeric vector of length \eqn{K-1}.
#' @keywords internal
.estimate_sigma2_alr <- function(R, X_std, B_minus, ref,
                                 delta = 1e-2,
                                 P_eff = ncol(X_std),
                                 floor = 1e-8) {
  N <- nrow(R); K <- ncol(R)
  nr <- setdiff(seq_len(K), ref)
  Rs <- R + delta
  Rs <- Rs / rowSums(Rs)
  ALR_R <- log(Rs[, nr, drop = FALSE] / Rs[, ref])
  mu_hat <- X_std %*% B_minus                    ## N x (K-1)
  resid <- ALR_R - mu_hat
  df <- max(N - P_eff, 1L)
  s2 <- colSums(resid^2) / df
  pmax(s2, floor)
}
