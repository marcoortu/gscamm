## ---------------------------------------------------------------------------
## Utility functions: standardization, perplexity, alignment, RMSE.
## ---------------------------------------------------------------------------

#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Column standardization (mean 0, sd 1) returning attributes for predict
#'
#' @param X numeric matrix \eqn{N \times P}.
#' @keywords internal
.standardize <- function(X) {
  mu <- colMeans(X)
  sd_ <- apply(X, 2, stats::sd)
  ## constant columns: keep them centered with sd 1 to avoid division by 0
  sd_[sd_ == 0] <- 1
  Xs <- sweep(X, 2, mu, "-")
  Xs <- sweep(Xs, 2, sd_, "/")
  attr(Xs, "center") <- mu
  attr(Xs, "scale")  <- sd_
  Xs
}

#' Perplexity of a fitted mixture model on a count matrix
#'
#' \eqn{\mathrm{Perplexity} = \exp\{ -\sum_{i,v} w_{iv} \log
#' (\Theta \Phi)_{iv} / \sum_{i,v} w_{iv} \}}.
#'
#' @param W observed count matrix \eqn{N \times V}.
#' @param Theta mixture proportions \eqn{N \times K}.
#' @param Phi component-category distributions \eqn{K \times V}.
#'
#' @return positive scalar.
#' @export
perplexity <- function(W, Theta, Phi) {
  q <- Theta %*% Phi
  q[q < .Machine$double.xmin] <- .Machine$double.xmin
  total <- sum(W)
  if (total == 0) return(NA_real_)
  exp(-sum(W * log(q)) / total)
}

#' Align estimated components to a reference via the Hungarian algorithm
#'
#' Given two component-by-category matrices, finds the permutation of
#' the rows of \code{Phi_hat} that minimizes the total \eqn{L_1} distance
#' to \code{Phi_ref}. Standard practice for label-switching in mixture
#' models, used in the simulation study of the paper.
#'
#' @param Phi_hat estimated component-category matrix \eqn{K \times V}.
#' @param Phi_ref reference component-category matrix \eqn{K \times V}.
#'
#' @return integer permutation \code{p} of length \eqn{K} such that
#'   \code{Phi_hat[p, ]} is aligned with \code{Phi_ref}.
#' @export
align_components <- function(Phi_hat, Phi_ref) {
  K1 <- nrow(Phi_hat); K2 <- nrow(Phi_ref)
  if (K1 != K2)
    stop("align_components requires matrices with the same number of rows.")
  D <- matrix(0, K1, K2)
  for (a in seq_len(K1))
    for (b in seq_len(K2))
      D[a, b] <- sum(abs(Phi_hat[a, ] - Phi_ref[b, ]))
  ## solve_LSAP minimizes sum of D[i, p[i]]; we want for each ref column b
  ## the row a of Phi_hat that should be relabeled to b.
  perm <- as.integer(clue::solve_LSAP(D))
  ## perm[a] = b means row a of Phi_hat -> position b in the aligned matrix.
  ## We want the inverse: which row of Phi_hat to place at position b.
  inv <- integer(K1)
  inv[perm] <- seq_len(K1)
  inv
}

#' RMSE between two matrices
#'
#' @param A,B matrices of the same shape.
#' @return non-negative scalar.
#' @keywords internal
.rmse <- function(A, B) sqrt(mean((A - B)^2))

#' Apply a permutation to the rows of Phi and the columns of Theta / B
#'
#' Reorders the components of a fitted GSCA-MM object so that
#' \code{Phi_new[k, ] = Phi_old[perm[k], ]}, with the matching column
#' permutations on \code{Theta} and \code{B}.
#'
#' @param fit a \code{gscamm} object.
#' @param perm integer vector, permutation of \code{1:K}.
#' @return a \code{gscamm} object with permuted components.
#' @keywords internal
.permute_fit <- function(fit, perm) {
  fit$Phi   <- fit$Phi[perm, , drop = FALSE]
  fit$Theta <- fit$Theta[, perm, drop = FALSE]
  fit$B     <- fit$B[, perm, drop = FALSE]
  ## Gamma component-score block (last K columns) must be permuted too
  P <- ncol(fit$X_std)
  K <- length(perm)
  fit$Gamma[, (P + 1):(P + K)] <- fit$Gamma[, P + perm, drop = FALSE]
  fit
}
