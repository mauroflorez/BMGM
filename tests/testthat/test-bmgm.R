test_that("bmgm runs without errors", {
  X <- cbind(rnorm(100), rpois(100, 5), rnorm(100, mean = 2), rbinom(100, 5, 0.8))
  type <- c("c", "d", "c", "d")
  expect_silent(BMGM::bmgm(X, type, nburn = 100, nsample = 100))
})

test_that("bmgm handles multinomial categorical variables", {
  set.seed(42)
  n <- 100
  X1 <- rpois(n, lambda = 3)
  X2 <- rnorm(n, mean = 5)
  X3 <- rnorm(n, mean = X2)
  X4 <- sample(x = c(1, 2, 3, 4), size = n, replace = TRUE, prob = rep(1/4, 4))
  X <- cbind(X1, X2, X3, X4)
  type <- c("d", "c", "c", "m")
  expect_no_error(BMGM::bmgm(X, type, nburn = 100, nsample = 100))
})
