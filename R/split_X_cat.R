#' Split Categorical Variables
#'
#' @param X Input matrix
#' @param type Vector of variable types
#' @return A list with design matrix and category sizes
#' @export
split_X_cat <- function(X, type) {
  p <- ncol(X)
  X_design <- do.call(cbind, lapply(1:p, function(s) {
    if (type[s] == "m") stats::model.matrix(~factor(X[,s]))[, -1] else X[,s]
  }))
  categories <- ifelse(type == "m", apply(X, 2, function(x) length(unique(x)) - 1), 1)

  return(list("matrix" = X_design, "categories" = categories))
}
