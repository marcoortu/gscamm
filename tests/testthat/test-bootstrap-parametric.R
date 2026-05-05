test_that("parametric bootstrap returns valid CIs and B_alr matches plug-in", {
  sim <- simulate_gscamm(N = 80, V = 50, K = 3, P = 2, seed = 21,
                         doc_length_mean = 300)
  fit <- fit_gscamm(sim$W, sim$X, K = 3,
                    control = gscamm_control(max_iter = 30, tol = 1e-3),
                    seed = 21)
  eff_pi <- covariate_effects(fit)
  eff_pb <- bootstrap_covariate_effects(fit, B = 30, method = "parametric",
                                        type = "basic", seed = 1)

  expect_s3_class(eff_pb, "gscamm_effects_boot")
  expect_equal(eff_pb$resample_method, "parametric")
  expect_equal(eff_pb$type, "basic")
  expect_equal(dim(eff_pb$B_alr), dim(eff_pi$B_alr))

  ## CI bounds finite and lo <= hi
  co <- eff_pb$coefficients
  expect_true(all(is.finite(co$conf.low)))
  expect_true(all(is.finite(co$conf.high)))
  expect_true(all(co$conf.low <= co$conf.high))

  ## B_alr point estimate is bit-identical to the plug-in (bootstrap is
  ## strictly additive, must not modify the original fit's coefficients)
  expect_equal(eff_pb$B_alr, eff_pi$B_alr, tolerance = 0)
})

test_that("default method is parametric and back-compat infers noise_augmented", {
  sim <- simulate_gscamm(N = 60, V = 40, K = 3, P = 2, seed = 33,
                         doc_length_mean = 250)
  fit <- fit_gscamm(sim$W, sim$X, K = 3,
                    control = gscamm_control(max_iter = 20, tol = 1e-3),
                    seed = 33)

  ## default: parametric
  eff_default <- bootstrap_covariate_effects(fit, B = 15, seed = 1)
  expect_equal(eff_default$resample_method, "parametric")

  ## back-compat: noise_augment=TRUE without method= -> noise_augmented
  eff_legacy <- bootstrap_covariate_effects(fit, B = 15,
                                            noise_augment = TRUE, seed = 1)
  expect_equal(eff_legacy$resample_method, "noise_augmented")

  ## explicit method= overrides noise_augment
  eff_explicit <- bootstrap_covariate_effects(fit, B = 15,
                                              method = "parametric",
                                              noise_augment = TRUE, seed = 1)
  expect_equal(eff_explicit$resample_method, "parametric")
})

test_that("parametric bootstrap CIs widen with smaller L (more multinomial noise)", {
  ## With shorter documents the multinomial noise dominates; CIs should be
  ## strictly wider than with very long documents (signal-rich regime).
  sim_short <- simulate_gscamm(N = 80, V = 40, K = 3, P = 2, seed = 41,
                               doc_length_mean = 50)
  sim_long  <- simulate_gscamm(N = 80, V = 40, K = 3, P = 2, seed = 41,
                               doc_length_mean = 1000)
  ctl <- gscamm_control(max_iter = 30, tol = 1e-3)
  f_s <- fit_gscamm(sim_short$W, sim_short$X, K = 3, control = ctl, seed = 41)
  f_l <- fit_gscamm(sim_long$W,  sim_long$X,  K = 3, control = ctl, seed = 41)

  e_s <- bootstrap_covariate_effects(f_s, B = 30, method = "parametric",
                                     type = "basic", seed = 7)
  e_l <- bootstrap_covariate_effects(f_l, B = 30, method = "parametric",
                                     type = "basic", seed = 7)

  width_s <- mean(e_s$coefficients$conf.high - e_s$coefficients$conf.low)
  width_l <- mean(e_l$coefficients$conf.high - e_l$coefficients$conf.low)
  expect_gt(width_s, width_l)
})

test_that("seed reproducibility: same seed -> identical bootstrap draws", {
  sim <- simulate_gscamm(N = 60, V = 30, K = 3, P = 2, seed = 51,
                         doc_length_mean = 200)
  fit <- fit_gscamm(sim$W, sim$X, K = 3,
                    control = gscamm_control(max_iter = 20, tol = 1e-3),
                    seed = 51)
  e1 <- bootstrap_covariate_effects(fit, B = 15, method = "parametric",
                                    type = "basic", seed = 99)
  e2 <- bootstrap_covariate_effects(fit, B = 15, method = "parametric",
                                    type = "basic", seed = 99)
  expect_equal(e1$B_draws, e2$B_draws)
})
