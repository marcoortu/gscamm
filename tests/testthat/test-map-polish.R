test_that("MAP polish returns valid simplex output", {
  sim <- simulate_gscamm(N = 80, V = 60, K = 4, P = 3, seed = 11,
                         doc_length_mean = 400)
  fit <- fit_gscamm(sim$W, sim$X, K = 4,
                    control = gscamm_control(max_iter = 30, tol = 1e-4),
                    seed = 11)
  expect_false(is.null(fit$Theta_map))
  expect_equal(dim(fit$Theta_map), c(80, 4))
  expect_equal(rowSums(fit$Theta_map), rep(1, 80), tolerance = 1e-10)
  expect_true(all(fit$Theta_map >= 0))
  expect_equal(length(fit$sigma2_polish), 3)
  expect_true(all(fit$sigma2_polish > 0))
})

test_that("polish='none' produces null Theta_map", {
  sim <- simulate_gscamm(N = 60, V = 40, K = 3, P = 2, seed = 17,
                         doc_length_mean = 300)
  fit <- fit_gscamm(sim$W, sim$X, K = 3, polish = "none",
                    control = gscamm_control(max_iter = 20, tol = 1e-3),
                    seed = 17)
  expect_null(fit$Theta_map)
  expect_null(fit$sigma2_polish)
  expect_false(is.null(fit$sigma2_k))   ## diagnostic still computed
})

test_that("MAP polish reduces structural-Theta error when structural is sub-optimal", {
  ## Recreate a scaled-down version of the paper simulation regime where
  ## the structural Theta is notably biased relative to the noisy truth
  ## (rmse_theta in the paper sims is ~0.15). We expect Theta_map to
  ## materially improve over Theta_struct here.
  sim <- simulate_gscamm(N = 400, V = 300, K = 10, P = 6, seed = 31,
                         doc_length_mean = 600, sigma = 0.3)
  fit <- fit_gscamm(sim$W, sim$X, K = 10,
                    control = gscamm_control(max_iter = 100, tol = 1e-5),
                    seed = 31)
  perm <- align_components(fit$Phi, sim$Phi)
  Theta_a     <- fit$Theta[, perm]
  Theta_map_a <- fit$Theta_map[, perm]
  rmse_struct <- sqrt(mean((Theta_a     - sim$Theta)^2))
  rmse_map    <- sqrt(mean((Theta_map_a - sim$Theta)^2))
  ## At the paper-sim regime, Theta_struct ~ 0.13 and the polish should
  ## bring it down by at least 30%.
  expect_gt(rmse_struct, 0.10)
  expect_lt(rmse_map, 0.7 * rmse_struct)
})

test_that("MAP polish degenerates to structural Theta as sigma2 -> 0", {
  sim <- simulate_gscamm(N = 60, V = 40, K = 3, P = 2, seed = 41,
                         doc_length_mean = 300)
  fit <- fit_gscamm(sim$W, sim$X, K = 3, sigma2_polish = 1e-6,
                    control = gscamm_control(max_iter = 20, tol = 1e-3),
                    seed = 41)
  expect_lt(max(abs(fit$Theta_map - fit$Theta)), 1e-2)
})

test_that("sigma2_polish accepts scalar, vector, and 'auto'", {
  sim <- simulate_gscamm(N = 60, V = 40, K = 3, P = 2, seed = 53,
                         doc_length_mean = 300)
  ctl <- gscamm_control(max_iter = 15, tol = 1e-3)

  fit_scalar <- fit_gscamm(sim$W, sim$X, K = 3,
                           sigma2_polish = 0.5, control = ctl, seed = 53)
  expect_equal(unname(fit_scalar$sigma2_polish), rep(0.5, 2))

  fit_vec <- fit_gscamm(sim$W, sim$X, K = 3,
                        sigma2_polish = c(0.1, 0.4),
                        control = ctl, seed = 53)
  expect_equal(unname(fit_vec$sigma2_polish), c(0.1, 0.4))

  fit_auto <- fit_gscamm(sim$W, sim$X, K = 3, sigma2_polish = "auto",
                         sigma2_max = 0.8,
                         control = ctl, seed = 53)
  expect_true(all(fit_auto$sigma2_polish <= 0.8 + 1e-12))
  expect_true(all(fit_auto$sigma2_polish > 0))
})

test_that("MAP polish does not change B or its sandwich CIs", {
  sim <- simulate_gscamm(N = 80, V = 50, K = 3, P = 2, seed = 67,
                         doc_length_mean = 300)
  ctl <- gscamm_control(max_iter = 30, tol = 1e-4)
  fit_with    <- fit_gscamm(sim$W, sim$X, K = 3, polish = "map",
                            control = ctl, seed = 67)
  fit_without <- fit_gscamm(sim$W, sim$X, K = 3, polish = "none",
                            control = ctl, seed = 67)
  expect_equal(fit_with$B, fit_without$B)
  eff_with    <- covariate_effects(fit_with)
  eff_without <- covariate_effects(fit_without)
  expect_equal(eff_with$B_alr, eff_without$B_alr)
  expect_equal(eff_with$coefficients$std.error,
               eff_without$coefficients$std.error,
               tolerance = 1e-10)
})
