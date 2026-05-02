test_that("covariate_effects returns expected structure", {
  sim <- simulate_gscamm(N = 80, V = 30, K = 4, P = 3, seed = 5,
                         doc_length_mean = 200)
  fit <- fit_gscamm(sim$W, sim$X, K = 4,
                    control = gscamm_control(max_iter = 25), seed = 5)
  eff <- covariate_effects(fit)
  expect_s3_class(eff, "gscamm_effects")
  expect_equal(dim(eff$B_alr), c(fit$P + 1, fit$K - 1))
  expect_true(all(c("estimate", "std.error", "p.value", "p.adj",
                    "odds.ratio", "or.low", "or.high") %in%
                  names(eff$coefficients)))
  ## SEs are non-negative
  expect_true(all(eff$coefficients$std.error >= 0, na.rm = TRUE))
  ## CIs are ordered
  expect_true(all(eff$coefficients$conf.low <= eff$coefficients$conf.high))
})

test_that("alignment finds the identity for an aligned matrix", {
  set.seed(0)
  P <- matrix(runif(15), 3, 5)
  P <- P / rowSums(P)
  perm <- align_components(P, P)
  expect_equal(perm, 1:3)
})
