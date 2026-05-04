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
  if (!is.null(fit$R)) fit$R <- fit$R[, perm, drop = FALSE]
  fit$B     <- fit$B[, perm, drop = FALSE]
  ## Gamma component-score block (last K columns) must be permuted too
  P <- ncol(fit$X_std)
  K <- length(perm)
  fit$Gamma[, (P + 1):(P + K)] <- fit$Gamma[, P + perm, drop = FALSE]
  ## refresh the canonical non-reference path-coefficient block. This is
  ## valid only when perm does not move the reference component, which is
  ## the contract of the alignment routine used by the bootstrap.
  ref <- fit$gsca_ref %||% K
  fit$B_minus <- fit$B[, -ref, drop = FALSE]
  rn <- rownames(fit$B_minus)
  if (is.null(rn) && !is.null(colnames(fit$X_std))) rn <- colnames(fit$X_std)
  if (is.null(rn)) rn <- paste0("X", seq_len(nrow(fit$B_minus)))
  dimnames(fit$B_minus) <- list(rn,
                                paste0("comp", setdiff(seq_len(K), ref)))
  fit
}

#' k-means++ seeded initialization for the component-category matrix Phi
#'
#' Builds a non-pathological starting Phi by k-means++ seeding the
#' row-normalized count matrix and refining the K centers with a few Lloyd
#' iterations on a subsample of rows. Reduces the local-optima frequency
#' of the EM-GSCA loop relative to a Dirichlet random start, especially
#' with K >= 5 components.
#'
#' Cost is O(K * N_sub * V) for the seeding pass plus a constant number
#' of Lloyd iterations on at most \code{subsample} rows. Designed to add
#' negligible overhead to \code{\link{fit_gscamm}} even when called inside
#' a bootstrap loop.
#'
#' @param W non-negative count matrix \eqn{N \times V}.
#' @param K number of components.
#' @param subsample maximum number of rows used for the kmeans iterations
#'   (default 500).
#' @param iter.max Lloyd iterations after seeding (default 10).
#' @return numeric matrix \eqn{K \times V} on the simplex.
#' @keywords internal
.init_phi_kmeans <- function(W, K, subsample = 500L, iter.max = 10L) {
  N <- nrow(W); V <- ncol(W)
  rs <- pmax(rowSums(W), 1)
  P <- W / rs                                       ## N x V row-stochastic

  ## restrict to a subsample for cost control on large N
  if (N > subsample) {
    idx_sub <- sample.int(N, subsample)
    Psub <- P[idx_sub, , drop = FALSE]
  } else {
    Psub <- P
  }
  Nsub <- nrow(Psub)

  ## degenerate edge case: fewer non-empty rows than K -> fall back to
  ## Dirichlet random init at the caller
  nz <- which(rowSums(Psub) > 0)
  if (length(nz) < K) {
    Phi <- matrix(stats::rgamma(K * V, 0.1, 1), K, V)
    return(Phi / rowSums(Phi))
  }

  ## ---- k-means++ seeding on Psub --------------------------------------
  seeds <- integer(K)
  seeds[1] <- sample(nz, 1L)
  ## squared L1 distances (cheap proxy for sparse counts)
  d <- rowSums(abs(sweep(Psub, 2, Psub[seeds[1], ], "-")))
  d2 <- d^2
  for (k in 2:K) {
    if (sum(d2) <= 0) {
      seeds[k] <- sample(setdiff(nz, seeds[seq_len(k - 1L)]), 1L)
    } else {
      seeds[k] <- sample.int(Nsub, 1L, prob = d2)
    }
    new_d <- rowSums(abs(sweep(Psub, 2, Psub[seeds[k], ], "-")))
    d <- pmin(d, new_d)
    d2 <- d^2
  }
  centers <- Psub[seeds, , drop = FALSE]            ## K x V

  ## ---- a few Lloyd iterations to refine centers -----------------------
  if (iter.max > 0L) {
    for (it in seq_len(iter.max)) {
      ## assign each row in Psub to nearest center (L1)
      ## D_{ik} = sum_v |Psub_iv - centers_kv|
      D <- matrix(0, Nsub, K)
      for (k in seq_len(K))
        D[, k] <- rowSums(abs(sweep(Psub, 2, centers[k, ], "-")))
      a <- max.col(-D, ties.method = "first")
      ## update centers as cluster means; empty clusters keep previous value
      new_centers <- centers
      for (k in seq_len(K)) {
        rows_k <- which(a == k)
        if (length(rows_k) > 0L)
          new_centers[k, ] <- colMeans(Psub[rows_k, , drop = FALSE])
      }
      delta <- max(abs(new_centers - centers))
      centers <- new_centers
      if (delta < 1e-6) break
    }
  }

  ## smooth zeros and normalize to the simplex
  Phi <- centers + 1e-6
  Phi / rowSums(Phi)
}
