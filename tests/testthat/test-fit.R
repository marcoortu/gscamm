test_that("fit_gscamm returns a well-formed gscamm object", {
  sim <- simulate_gscamm(N = 60, V = 30, K = 3, P = 3, seed = 42,
                         doc_length_mean = 200)
  fit <- fit_gscamm(sim$W, sim$X, K = 3,
                    control = gscamm_control(max_iter = 20, tol = 1e-3),
                    seed = 1)
  expect_s3_class(fit, "gscamm")
  expect_equal(dim(fit$Phi), c(3, 30))
  expect_equal(dim(fit$Theta), c(60, 3))
  expect_equal(dim(fit$B), c(3, 3))
  ## simplex constraints
  expect_equal(rowSums(fit$Theta), rep(1, 60))
  expect_equal(rowSums(fit$Phi),   rep(1, 3))
  ## perplexity is finite
  expect_true(is.finite(tail(fit$convergence$perplexity, 1)))
})

test_that("perplexity trajectory stabilizes across iterations", {
  ## EM-GSCA is not strictly perplexity-monotone (the GSCA ridge step is
  ## not a perplexity M-step), so we test that the trajectory stabilizes
  ## at a value close to the run minimum rather than strict monotonicity.
  sim <- simulate_gscamm(N = 80, V = 50, K = 4, P = 3, seed = 7,
                         doc_length_mean = 150)
  fit <- fit_gscamm(sim$W, sim$X, K = 4,
                    control = gscamm_control(max_iter = 30, tol = 1e-5,
                                             trace_perplexity = TRUE),
                    seed = 7)
  perp <- fit$convergence$perplexity
  ## last 5 iterations are close to each other (relative range < 1%)
  last5 <- tail(perp, 5)
  expect_lt(diff(range(last5)) / mean(last5), 0.01)
  ## final value is within 1% of the run minimum
  expect_lt((perp[length(perp)] - min(perp)) / min(perp), 0.01)
})

test_that("predict returns valid simplex matrices", {
  sim <- simulate_gscamm(N = 50, V = 25, K = 3, P = 2, seed = 11,
                         doc_length_mean = 200)
  fit <- fit_gscamm(sim$W, sim$X, K = 3,
                    control = gscamm_control(max_iter = 15), seed = 11)
  Th <- predict(fit, sim$X[1:5, , drop = FALSE])
  expect_equal(dim(Th), c(5, 3))
  expect_equal(rowSums(Th), rep(1, 5))
})

test_that("link choices all run end-to-end", {
  sim <- simulate_gscamm(N = 40, V = 20, K = 3, P = 2, seed = 99,
                         doc_length_mean = 100)
  for (lk in c("logistic_normal", "dirichlet", "zero_inflated")) {
    fit <- fit_gscamm(sim$W, sim$X, K = 3, link = lk,
                      control = gscamm_control(max_iter = 10), seed = 1)
    expect_s3_class(fit, "gscamm")
    expect_equal(rowSums(fit$Theta), rep(1, 40))
  }
})

test_that("identifiability constraint: B[, ref] is structurally zero", {
  sim <- simulate_gscamm(N = 60, V = 25, K = 4, P = 3, seed = 13,
                         doc_length_mean = 200)
  for (sp in c("alr", "simplex")) {
    fit <- fit_gscamm(sim$W, sim$X, K = 4, gsca_space = sp,
                      control = gscamm_control(max_iter = 15), seed = 13)
    expect_equal(max(abs(fit$B[, fit$gsca_ref])), 0,
                 info = paste("ref column not zero for gsca_space=", sp))
    expect_equal(dim(fit$B_minus), c(3L, 3L))
  }
})

test_that("coef.gscamm returns the canonical block by default", {
  sim <- simulate_gscamm(N = 50, V = 20, K = 4, P = 3, seed = 21,
                         doc_length_mean = 200)
  fit <- fit_gscamm(sim$W, sim$X, K = 4,
                    control = gscamm_control(max_iter = 12), seed = 21)
  expect_equal(dim(coef(fit)),                c(3L, 3L))
  expect_equal(dim(coef(fit, augmented = TRUE)), c(3L, 4L))
})

test_that("default gsca_space is alr", {
  sim <- simulate_gscamm(N = 30, V = 15, K = 3, P = 2, seed = 99,
                         doc_length_mean = 100)
  fit <- fit_gscamm(sim$W, sim$X, K = 3,
                    control = gscamm_control(max_iter = 5), seed = 1)
  expect_identical(fit$gsca_space, "alr")
})
