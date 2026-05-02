test_that("link functions return rows on the simplex", {
  S <- matrix(rnorm(6 * 4), 6, 4)
  for (g in list(link_logistic_normal, link_dirichlet,
                 function(x) link_zero_inflated(x, eps = 1e-3))) {
    Th <- g(S)
    expect_equal(rowSums(Th), rep(1, nrow(S)))
    expect_true(all(Th >= 0))
  }
})

test_that("zero_inflated produces some exact zeros for low scores", {
  set.seed(1)
  S <- matrix(rnorm(20 * 5, sd = 3), 20, 5)
  Th <- link_zero_inflated(S, eps = 0.1)
  expect_true(any(Th == 0))
  expect_equal(rowSums(Th), rep(1, nrow(S)))
})
