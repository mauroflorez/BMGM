#' Split Categorical Variables
#'
#' @param X Input matrix
#' @param type Vector of variable types
#' @return A list with design matrix and category sizes
#' @export
split_X_cat <- function(X, type, categories = NULL) {
  p <- ncol(X)

  # Compute categories from data if not provided
  if(is.null(categories)){
    categories <- ifelse(type == "m",
                         apply(X, 2, function(x) length(unique(x[!is.na(x)])) - 1),
                         1)
  }

  X_design <- do.call(cbind, lapply(1:p, function(s) {
    if (type[s] == "m") {
      # Use fixed levels to ensure consistent dummy columns even if some
      # categories are absent after imputation
      K <- categories[s] + 1  # total number of categories (including reference)
      stats::model.matrix(~factor(X[,s], levels = 1:K))[, -1, drop = FALSE]
    } else {
      X[,s]
    }
  }))

  return(list("matrix" = X_design, "categories" = categories))
}
