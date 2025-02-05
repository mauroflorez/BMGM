test_that("bmgm runs without errors", {
  X <- cbind(rnorm(100), rpois(100, 5), rnorm(100, mean = 2), rbinom(100, 5, 0.8))
  type <- c("c", "d", "c", "d")
  expect_silent(BMGM::bmgm(X, type, nburn = 100, nsample = 100))
})
