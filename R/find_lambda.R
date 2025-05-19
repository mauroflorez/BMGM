#' Find Optimal Lambda for Transformation
#'
#' @param X Input data
#' @param type Vector of variable types
#' @return Optimal lambda estimate
#' @export
find_lambda <- function(X, type) {
  d <- which(type != "m")
  X_d <- X[, d]
  p0 <- max(X_d, na.rm = T)
  p_v <- seq(0, 2 * p0, by = 0.1)

  errors <- sapply(p_v, function(i) {
    F_X <- F_transformation(X_d, type[d], parameter = i, cont = FALSE)
    sqrt(sum((cov(X_d) - cov(F_X))^2))
  })

  return(p_v[which.min(errors)])
}
