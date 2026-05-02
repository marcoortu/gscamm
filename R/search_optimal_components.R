## ---------------------------------------------------------------------------
## Component selection: fits the model for each K in a grid and reports
## perplexity (and optionally a held-out variant). Generalizes the
## search_optimal_topics utility mentioned in Section 5 of the paper.
## ---------------------------------------------------------------------------

#' Search the optimal number of components by perplexity
#'
#' Fits the GSCA-MM model for each \eqn{K} in a user-provided grid and
#' reports perplexity (in-sample, and optionally held-out via random
#' sub-sampling of the count matrix). The user can then pick the \eqn{K}
#' that minimizes (or stabilizes) the chosen criterion.
#'
#' @param W count matrix.
#' @param X covariate matrix or data frame.
#' @param Ks integer vector of candidate component counts.
#' @param link see \code{\link{fit_gscamm}}.
#' @param control see \code{\link{gscamm_control}}.
#' @param holdout fraction of token mass held out per row for held-out
#'   perplexity (default 0, meaning in-sample only).
#' @param seed optional integer for reproducibility (governs the holdout
#'   sampling and the model initialization).
#' @param verbose logical.
#'
#' @return data frame with columns \code{K}, \code{perplexity}, and
#'   optionally \code{holdout_perplexity}.
#' @export
search_optimal_components <- function(W, X, Ks,
                                      link = "logistic_normal",
                                      control = gscamm_control(),
                                      holdout = 0,
                                      seed = NULL,
                                      verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  W <- as.matrix(W)
  ## hold-out split: for each cell, keep a Binomial(n_iv, 1-holdout) draw
  if (holdout > 0) {
    if (holdout >= 1) stop("`holdout` must be in [0, 1).")
    Wh <- matrix(stats::rbinom(length(W), as.integer(W), holdout),
                 nrow(W), ncol(W))
    Wt <- W - Wh
  } else {
    Wt <- W; Wh <- NULL
  }

  out <- data.frame(K = integer(0), perplexity = double(0))
  if (!is.null(Wh)) out$holdout_perplexity <- double(0)

  for (K in Ks) {
    if (verbose) message("Fitting K = ", K)
    fit <- fit_gscamm(Wt, X, K = K, link = link, control = control,
                      verbose = FALSE, seed = seed)
    perp_in <- perplexity(Wt, fit$Theta, fit$Phi)
    row <- data.frame(K = K, perplexity = perp_in)
    if (!is.null(Wh)) {
      row$holdout_perplexity <- perplexity(Wh, fit$Theta, fit$Phi)
    }
    out <- rbind(out, row)
  }
  out
}
