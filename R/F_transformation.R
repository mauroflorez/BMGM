#' Transformation Function for Mixed Graphical Models
#'
#' @description Performs transformation of input data based on its type.
#' Used for continuous, discrete, and zero-inflated data in Mixed Graphical Models.
#'
#' @param X A numeric matrix or vector.
#' @param type A character vector specifying the type of each variable.
#'        Options: "c" (continuous), "d" (discrete), "z" (zero-inflated), "m" (categorical).
#' @param parameter Numeric parameter used in the transformation.
#' @param cont Logical. If TRUE, continuous variables are transformed. Default is FALSE.
#'
#' @return A matrix of transformed data.
#' @export
F_transformation <- function(X, type, parameter, cont = FALSE){
  X_t <- X
  transform <- function(x, typ, param){
    switch(typ,
           m = x,
           d = atan(x)^param,
           z = atan(x)^param,
           c = if(cont) sign(x)*atan(abs(x))^param else x)
  }

  if(!(is.matrix(X) | is.data.frame(X))) X_t <- transform(X, type, parameter) else
    for(c in 1:ncol(X)) X_t[,c] <- transform(X[,c], type[c], parameter)

  return(X_t)
}
