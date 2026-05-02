## ---------------------------------------------------------------------------
## Main fit function for the GSCA-MM framework. Implements Algorithm 1 of
## Ortu and Frigau (2026), generalized from topic models to arbitrary
## mixture-of-multinomials models with unit-level covariates.
## ---------------------------------------------------------------------------

#' Control parameters for the EM-GSCA algorithm
#'
#' @param max_iter maximum number of EM-GSCA iterations (default 100).
#' @param tol convergence tolerance on the change in
#'   \eqn{(\boldsymbol{\Phi}, \boldsymbol{\Theta}, \boldsymbol{\Gamma})}
#'   measured as the maximum absolute change across the three matrices
#'   (default \code{1e-4}).
#' @param alpha Dirichlet smoothing constant added to the sufficient
#'   statistics in the M-step (default \code{0.01}).
#' @param lambda0 initial ridge parameter \eqn{\lambda_0} for the geometric
#'   schedule \eqn{\lambda_B^{(t)} = \lambda_0 \rho^t} (default 1).
#' @param rho ridge decay rate \eqn{\rho \in (0, 1)} (default 0.9).
#' @param eps sparsity threshold used by the zero-inflated link
#'   (default \code{1e-3}).
#' @param init_Phi optional list with \code{Phi} (\eqn{K \times V}) for
#'   warm-starting; default initializes from a Dirichlet on the simplex.
#' @param init_B optional warm-start \eqn{P \times K} matrix; default zero.
#' @param min_iter minimum number of iterations before convergence checks
#'   are applied (default 5).
#' @param trace_perplexity if \code{TRUE}, perplexity is recomputed every
#'   iteration and stored; otherwise only every \code{trace_every}
#'   iterations.
#' @param trace_every interval at which to record perplexity when
#'   \code{trace_perplexity = FALSE} (default 5).
#'
#' @return a list of class \code{gscamm_control}.
#' @export
gscamm_control <- function(max_iter = 100,
                           tol = 1e-4,
                           alpha = 0.01,
                           lambda0 = 1,
                           rho = 0.9,
                           eps = 1e-3,
                           init_Phi = NULL,
                           init_B = NULL,
                           min_iter = 5,
                           trace_perplexity = TRUE,
                           trace_every = 1) {
  stopifnot(max_iter >= 1, tol >= 0, alpha >= 0,
            lambda0 >= 0, rho > 0, rho < 1, eps > 0,
            min_iter >= 1, trace_every >= 1)
  structure(list(max_iter = max_iter, tol = tol, alpha = alpha,
                 lambda0 = lambda0, rho = rho, eps = eps,
                 init_Phi = init_Phi, init_B = init_B,
                 min_iter = min_iter,
                 trace_perplexity = trace_perplexity,
                 trace_every = trace_every),
            class = "gscamm_control")
}

#' Fit a GSCA-MM model
#'
#' Fits the Generalized Structured Component Analysis for
#' Mixture-of-Multinomials (GSCA-MM) model of Ortu and Frigau (2026),
#' implementing the blockwise EM-GSCA algorithm (Algorithm 1) on a count
#' matrix \code{W} and a covariate matrix \code{X}.
#'
#' The model represents the count matrix as a mixture
#' \eqn{\mathbf{W} \approx \boldsymbol{\Theta}\boldsymbol{\Phi}}, with
#' mixture weights driven by linear combinations of the standardized
#' covariates through a path coefficient matrix \eqn{\mathbf{B}}. See
#' Section 2 of the paper for the full model specification.
#'
#' @param W non-negative numeric matrix \eqn{N \times V} of counts
#'   (e.g. document-term matrix, OTU table, basket-by-product matrix).
#' @param X numeric matrix or data frame \eqn{N \times P} of covariates.
#'   Standardized internally.
#' @param K positive integer, number of mixture components.
#' @param link link function mapping component scores to mixture
#'   proportions; one of \code{"logistic_normal"} (default),
#'   \code{"dirichlet"}, \code{"zero_inflated"}.
#' @param gsca_space geometry in which the GSCA step regresses the
#'   responsibilities on the covariates. \code{"simplex"} (default,
#'   matches Equation (8) of the paper) regresses raw responsibilities
#'   linearly on \code{X}; \code{"alr"} regresses
#'   \eqn{\log(r_{ik} / r_{i\,\mathrm{ref}})} on \code{X} with a
#'   user-chosen reference component \code{gsca_ref} (default \code{K}).
#'   The ALR variant aligns the linear regression with the natural
#'   geometry of compositional data and substantially reduces the
#'   structural bias of the path-coefficient estimator, at the cost of
#'   departing slightly from the original GSCA-MM specification.
#' @param gsca_ref integer reference component used when
#'   \code{gsca_space = "alr"} (default \code{K}).
#' @param control list of control parameters; see \code{\link{gscamm_control}}.
#' @param verbose logical, print per-iteration progress.
#' @param seed optional integer for reproducible initialization.
#'
#' @return an object of S3 class \code{gscamm} with components:
#'   \item{Phi}{\eqn{K \times V} component-category matrix.}
#'   \item{Theta}{\eqn{N \times K} mixture-proportion matrix.}
#'   \item{B}{\eqn{P \times K} path coefficient matrix.}
#'   \item{Gamma}{\eqn{N \times (P+K)} component-score matrix.}
#'   \item{X_std}{standardized covariate matrix used for fitting.}
#'   \item{W}{the input count matrix (kept for diagnostics).}
#'   \item{link, K, P, V, N}{model dimensions and link.}
#'   \item{convergence}{list with iteration history of perplexity, max
#'     change, and ridge schedule, plus convergence flag.}
#'   \item{control}{the resolved control list.}
#'
#' @references
#' Ortu, M. and Frigau, L. (2026). Structured Component Regression for
#' Covariate Effects in Mixture-of-Multinomials Models.
#'
#' @examples
#' \dontrun{
#'   sim <- simulate_gscamm(N = 200, V = 100, K = 5, P = 4)
#'   fit <- fit_gscamm(sim$W, sim$X, K = 5, link = "logistic_normal")
#'   eff <- covariate_effects(fit)
#'   summary(fit)
#' }
#' @export
fit_gscamm <- function(W, X, K,
                       link = c("logistic_normal", "dirichlet",
                                "zero_inflated"),
                       gsca_space = c("simplex", "alr"),
                       gsca_ref = K,
                       control = gscamm_control(),
                       verbose = FALSE,
                       seed = NULL) {
  link <- match.arg(link)
  gsca_space <- match.arg(gsca_space)
  if (gsca_ref < 1 || gsca_ref > K) stop("gsca_ref out of range.")
  if (!inherits(control, "gscamm_control"))
    control <- do.call(gscamm_control, as.list(control))
  if (!is.null(seed)) set.seed(seed)

  ## --- coerce inputs -------------------------------------------------------
  W <- as.matrix(W)
  if (!is.numeric(W) || any(W < 0))
    stop("W must be a non-negative numeric matrix.")
  if (is.data.frame(X)) X <- stats::model.matrix(~ . - 1, data = X)
  X <- as.matrix(X)
  if (nrow(W) != nrow(X))
    stop("W and X must have the same number of rows.")
  N <- nrow(W); V <- ncol(W); P <- ncol(X)
  if (K < 2) stop("K must be at least 2.")

  ## --- standardize covariates ----------------------------------------------
  X_std <- .standardize(X)
  XtX <- crossprod(X_std)

  ## --- initialize Phi on the simplex ---------------------------------------
  if (!is.null(control$init_Phi)) {
    Phi <- as.matrix(control$init_Phi)
    if (!identical(dim(Phi), c(K, V)))
      stop("init_Phi must have dimension K x V.")
  } else {
    ## random Dirichlet(0.1) rows
    Phi <- matrix(stats::rgamma(K * V, 0.1, 1), K, V)
    Phi <- Phi / rowSums(Phi)
  }

  ## --- initialize B and Gamma ----------------------------------------------
  if (!is.null(control$init_B)) {
    B <- as.matrix(control$init_B)
    if (!identical(dim(B), c(P, K)))
      stop("init_B must have dimension P x K.")
  } else {
    B <- matrix(0, P, K)
  }
  ## small random topic scores to break symmetry
  scores0 <- X_std %*% B + matrix(stats::rnorm(N * K, sd = 0.1), N, K)
  Theta <- .apply_link(scores0, link, eps = control$eps)
  Gamma <- cbind(X_std, scores0)

  ## --- iterate -------------------------------------------------------------
  hist_perp <- numeric(0)
  hist_dchange <- numeric(0)
  hist_lambda <- numeric(0)
  converged <- FALSE
  iter_done <- 0L

  for (t in seq_len(control$max_iter)) {
    iter_done <- t

    ## E-step
    R <- .e_step(W, Theta, Phi)

    ## M-step (Phi)
    Phi_new <- .m_step_phi(W, R, alpha = control$alpha)

    ## GSCA step (B): simplex- or ALR-space ridge
    lambda_t <- .lambda_schedule(t - 1L, control$lambda0, control$rho)
    B_new <- if (gsca_space == "simplex")
      .gsca_update_B(X_std, R, lambda_B = lambda_t, XtX = XtX)
    else
      .gsca_update_B_alr(X_std, R, lambda_B = lambda_t,
                         ref = gsca_ref, XtX = XtX)

    ## update Gamma and Theta from new scores
    scores_new <- X_std %*% B_new
    Theta_new <- .apply_link(scores_new, link, eps = control$eps)
    Gamma_new <- cbind(X_std, scores_new)

    ## convergence diagnostic: max abs change across (Phi, Theta, Gamma)
    dchange <- max(
      max(abs(Phi_new   - Phi)),
      max(abs(Theta_new - Theta)),
      max(abs(Gamma_new - Gamma))
    )

    ## perplexity tracking
    do_perp <- control$trace_perplexity ||
      (t %% control$trace_every == 0L) || (t == control$max_iter)
    if (do_perp) {
      perp <- perplexity(W, Theta_new, Phi_new)
      hist_perp <- c(hist_perp, perp)
    }
    hist_dchange <- c(hist_dchange, dchange)
    hist_lambda  <- c(hist_lambda, lambda_t)

    if (verbose) {
      msg <- sprintf("iter %3d  lambda_B = %.4g  d = %.4g", t, lambda_t, dchange)
      if (do_perp) msg <- paste0(msg, sprintf("  perplexity = %.3f", perp))
      message(msg)
    }

    ## commit
    Phi <- Phi_new; Theta <- Theta_new; B <- B_new; Gamma <- Gamma_new

    if (t >= control$min_iter && dchange < control$tol) {
      converged <- TRUE
      break
    }
  }

  ## final responsibilities R from the converged (Phi, Theta).  These
  ## reflect the data-driven posterior at the last E-step and are kept
  ## for the second-stage ALR-WLS regression (see covariate_effects()).
  R_final <- .e_step(W, Theta, Phi)

  fit <- list(
    Phi = Phi, Theta = Theta, R = R_final, B = B, Gamma = Gamma,
    X_std = X_std, X = X, W = W,
    link = link, gsca_space = gsca_space, gsca_ref = gsca_ref,
    K = K, P = P, V = V, N = N,
    convergence = list(
      iterations = iter_done,
      converged = converged,
      perplexity = hist_perp,
      d_change = hist_dchange,
      lambda_B = hist_lambda
    ),
    control = control,
    call = match.call()
  )
  class(fit) <- "gscamm"
  fit
}
