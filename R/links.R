#' Link functions mapping component scores to mixture proportions
#'
#' These functions implement the three link specifications discussed in
#' Section 2.1 of Ortu and Frigau (2026): logistic-normal, Dirichlet, and
#' zero-inflated. Each function takes the matrix of component scores
#' \eqn{\boldsymbol{\Gamma}^{(C)}} (the last \eqn{K} columns of the
#' component-score matrix \eqn{\boldsymbol{\Gamma}}) and returns the
#' mixture proportions \eqn{\boldsymbol{\Theta}}.
#'
#' @param scores numeric matrix of component scores, \eqn{N \times K}.
#' @param eps positive sparsity threshold for the zero-inflated link
#'   (default \code{1e-3}).
#'
#' @return numeric matrix \eqn{N \times K} of mixture proportions, with
#'   non-negative entries summing to one in each row.
#'
#' @name link_functions
NULL

#' @describeIn link_functions Logistic-normal (softmax) link, the default
#'   specification used in the paper. Equation (3).
#' @export
link_logistic_normal <- function(scores) {
  if (!is.matrix(scores)) scores <- as.matrix(scores)
  m <- apply(scores, 1, max)
  z <- exp(scores - m)
  z / rowSums(z)
}

#' @describeIn link_functions Dirichlet link with concentration parameters
#'   given by \eqn{\exp(\gamma_{ik})}. Equation (4). The deterministic
#'   point estimate used in the algorithm is the normalized mean of the
#'   Dirichlet, i.e. \code{exp(scores)} normalized rowwise; this matches
#'   the GSCA-MM estimation routine, which works with the conditional mean
#'   of the component scores given the covariates.
#' @export
link_dirichlet <- function(scores) {
  if (!is.matrix(scores)) scores <- as.matrix(scores)
  m <- apply(scores, 1, max)
  z <- exp(scores - m)
  z / rowSums(z)
}

#' @describeIn link_functions Zero-inflated link. Sets \eqn{\theta_{ik}=0}
#'   when \eqn{\exp(\gamma_{ik}) < \epsilon} and renormalizes the
#'   surviving entries. Equation (5).
#' @export
link_zero_inflated <- function(scores, eps = 1e-3) {
  if (!is.matrix(scores)) scores <- as.matrix(scores)
  m <- apply(scores, 1, max)
  z <- exp(scores - m)
  ## scale to "exp(gamma)" magnitude so that eps is comparable across rows;
  ## the threshold is applied to the *raw* exponentiated scores rescaled by
  ## the per-row max, which makes eps a relative cut-off.
  z[z < eps] <- 0
  rs <- rowSums(z)
  ## guard against rows with everything zeroed: fall back to plain softmax
  empty <- rs == 0
  if (any(empty)) {
    fb <- exp(scores[empty, , drop = FALSE] -
              apply(scores[empty, , drop = FALSE], 1, max))
    z[empty, ] <- fb
    rs[empty] <- rowSums(fb)
  }
  z / rs
}

## internal dispatcher
.apply_link <- function(scores, link, eps = 1e-3) {
  switch(link,
         logistic_normal = link_logistic_normal(scores),
         dirichlet       = link_dirichlet(scores),
         zero_inflated   = link_zero_inflated(scores, eps = eps),
         stop("Unknown link: ", link))
}
