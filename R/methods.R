## ---------------------------------------------------------------------------
## S3 methods for gscamm and gscamm_effects objects.
## ---------------------------------------------------------------------------

#' @export
print.gscamm <- function(x, ...) {
  cat("GSCA-MM fit\n")
  cat(sprintf("  observations N = %d, categories V = %d, components K = %d, covariates P = %d\n",
              x$N, x$V, x$K, x$P))
  cat(sprintf("  link        : %s\n", x$link))
  cat(sprintf("  iterations  : %d (%s)\n",
              x$convergence$iterations,
              if (x$convergence$converged) "converged" else "max_iter reached"))
  if (length(x$convergence$perplexity)) {
    cat(sprintf("  final perplexity : %.4f\n",
                x$convergence$perplexity[length(x$convergence$perplexity)]))
  }
  invisible(x)
}

#' @export
summary.gscamm <- function(object, top_n = 10, ...) {
  ## top categories per component
  V_names <- colnames(object$W)
  if (is.null(V_names)) V_names <- paste0("v", seq_len(object$V))
  top <- vector("list", object$K)
  for (k in seq_len(object$K)) {
    ord <- order(object$Phi[k, ], decreasing = TRUE)[seq_len(min(top_n, object$V))]
    top[[k]] <- data.frame(rank = seq_along(ord),
                           category = V_names[ord],
                           phi = object$Phi[k, ord],
                           stringsAsFactors = FALSE)
  }
  names(top) <- paste0("comp", seq_len(object$K))

  comp_share <- colMeans(object$Theta)
  out <- list(
    fit = object,
    top_categories = top,
    component_share = comp_share,
    perplexity = if (length(object$convergence$perplexity))
      object$convergence$perplexity[length(object$convergence$perplexity)]
      else NA_real_
  )
  class(out) <- "summary.gscamm"
  out
}

#' @export
print.summary.gscamm <- function(x, ...) {
  print(x$fit)
  cat("\nMean component share (colMeans of Theta):\n")
  print(round(x$component_share, 4))
  cat("\nTop categories per component:\n")
  for (k in seq_along(x$top_categories)) {
    cat(sprintf("\n  %s\n", names(x$top_categories)[k]))
    print(x$top_categories[[k]], row.names = FALSE)
  }
  invisible(x)
}

#' @export
coef.gscamm <- function(object, ...) {
  ## Returns the path coefficients on the standardized covariate scale.
  rn <- colnames(object$X_std); if (is.null(rn)) rn <- paste0("X", seq_len(object$P))
  cn <- paste0("comp", seq_len(object$K))
  dimnames(object$B) <- list(rn, cn)
  object$B
}

#' Predict mixture proportions for new observations
#'
#' Given a fitted GSCA-MM model and a new covariate matrix
#' \code{newdata}, computes predicted mixture proportions
#' \eqn{\hat{\theta}_i = g(\hat{B}^\top \mathbf{x}_i^{\mathrm{std}})}
#' under the same link \eqn{g} used at fitting. Standardization uses the
#' centering and scaling derived at fit time.
#'
#' @param object a \code{gscamm} fit.
#' @param newdata covariate matrix or data frame for new observations.
#' @param type one of \code{"theta"} (default, mixture proportions) or
#'   \code{"scores"} (component scores before applying the link).
#' @param ... unused.
#' @return numeric matrix.
#' @export
predict.gscamm <- function(object, newdata, type = c("theta", "scores"), ...) {
  type <- match.arg(type)
  if (is.data.frame(newdata)) newdata <- stats::model.matrix(~ . - 1, data = newdata)
  newdata <- as.matrix(newdata)
  if (ncol(newdata) != object$P)
    stop("newdata must have ", object$P, " columns to match the fitted model.")
  mu <- attr(object$X_std, "center")
  sc <- attr(object$X_std, "scale")
  Xs <- sweep(newdata, 2, mu, "-")
  Xs <- sweep(Xs, 2, sc, "/")
  scores <- Xs %*% object$B
  if (type == "scores") return(scores)
  .apply_link(scores, object$link, eps = object$control$eps)
}

#' Default plot of GSCA-MM diagnostics
#'
#' Plots perplexity and the maximum parameter change per iteration.
#'
#' @param x a \code{gscamm} fit.
#' @param ... passed to \code{plot}.
#' @export
plot.gscamm <- function(x, ...) {
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op), add = TRUE)
  perp <- x$convergence$perplexity
  if (length(perp))
    plot(seq_along(perp), perp, type = "b", pch = 19,
         xlab = "iteration (recorded)", ylab = "perplexity",
         main = "Perplexity")
  d <- x$convergence$d_change
  plot(seq_along(d), d, type = "b", pch = 19, log = "y",
       xlab = "iteration", ylab = "max abs change",
       main = "Convergence diagnostic")
  invisible(x)
}

#' @export
print.gscamm_effects <- function(x, digits = 4, ...) {
  cat("GSCA-MM covariate effects (ALR-WLS)\n")
  cat(sprintf("  reference component: %d\n", x$ref))
  cat(sprintf("  confidence level   : %.2f\n", x$level))
  cat(sprintf("  p-value adjustment : %s\n\n", x$adjust))
  print(format(x$coefficients, digits = digits))
  invisible(x)
}

#' Plot ALR coefficients with confidence intervals
#'
#' @param x a \code{gscamm_effects} object.
#' @param drop_intercept logical, drop intercept rows (default TRUE).
#' @param sig_only logical, restrict to significant adjusted p-values
#'   below 0.05 (default FALSE).
#' @param ... unused.
#' @export
plot.gscamm_effects <- function(x, drop_intercept = TRUE,
                                sig_only = FALSE, ...) {
  d <- x$coefficients
  if (drop_intercept) d <- d[d$covariate != "(Intercept)", , drop = FALSE]
  if (sig_only) d <- d[!is.na(d$p.adj) & d$p.adj < 0.05, , drop = FALSE]
  if (!nrow(d)) {
    message("No coefficients to plot.")
    return(invisible(x))
  }

  comps <- unique(d$component)
  op <- graphics::par(mfrow = c(1, length(comps)),
                      mar = c(4, 6, 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  for (cc in comps) {
    dd <- d[d$component == cc, , drop = FALSE]
    yl <- seq_len(nrow(dd))
    graphics::plot(dd$estimate, yl, xlim = range(c(dd$conf.low, dd$conf.high)),
                   yaxt = "n", xlab = "log-odds coefficient", ylab = "",
                   main = cc, pch = 19)
    graphics::axis(2, at = yl, labels = dd$covariate, las = 1)
    graphics::abline(v = 0, lty = 2, col = "grey50")
    graphics::segments(dd$conf.low, yl, dd$conf.high, yl)
  }
  invisible(x)
}
