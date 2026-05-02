#' gscamm: Generalized Structured Component Analysis for Mixture Models
#'
#' @description
#' Implementation of the GSCA-MM framework introduced by Ortu and Frigau
#' (2026), generalizing the GSCA-TM topic model to arbitrary
#' mixture-of-multinomials applications such as microbiome composition
#' analysis, market basket studies and topic modeling. Covariate effects
#' on mixture weights are represented through a path coefficient matrix
#' within a generalized structured component analysis formulation.
#' Estimation alternates EM updates for the component distributions with
#' a ridge-penalized regression step for the path coefficients,
#' decoupling latent structure recovery from structural parameter
#' estimation. Post-estimation inference on covariate effects is
#' available via additive log-ratio weighted least squares.
#'
#' @section Main entry points:
#' \describe{
#'   \item{\code{\link{fit_gscamm}}}{Fits the EM-GSCA algorithm
#'     (Algorithm 1 of the paper).}
#'   \item{\code{\link{covariate_effects}}}{Two-stage ALR-WLS inference
#'     on covariate effects (Section 2.4).}
#'   \item{\code{\link{simulate_gscamm}}}{Synthetic data generator
#'     mirroring the simulation design of Section 4.}
#'   \item{\code{\link{search_optimal_components}}}{Component-count
#'     selection by perplexity.}
#' }
#'
#' @keywords internal
"_PACKAGE"
