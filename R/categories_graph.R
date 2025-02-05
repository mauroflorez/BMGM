#' Convert Category Graph to General Graph
#'
#' @description This function converts a graph estimated from categorical variables
#' into a general adjacency matrix representation.
#'
#' @param q Number of transformed variables
#' @param p Number of original variables
#' @param var_names Vector of original variable names
#' @param esti_Z Adjacency matrix from context-specific inference
#' @param esti_Beta Estimated beta coefficients for graph edges
#' @param categories Vector indicating the number of categories for each variable
#'
#' @return A list with:
#' \item{Adj_Beta}{Estimated adjacency matrix for Beta}
#' \item{Adj_Z}{Estimated adjacency matrix for Z}
#' @export
categories_graph <- function(q, p, var_names, esti_Z, esti_Beta, categories){
  esti_Z_gen <- matrix(nrow = q, ncol = p)
  diag(esti_Z_gen) <- 0

  for(s in 1:p){
    if(categories[s] == 1){
      esti_Z_gen[,s] <- esti_Z[,which(var_names == s)]
    } else {
      esti_Z_gen[,s] <- pmin(1,rowSums(esti_Z[,which(var_names == s)]))
    }
  }

  esti_Z_gen_mod <- matrix(nrow = p, ncol = p)

  for(l in 1:p){
    if(categories[l] == 1){
      esti_Z_gen_mod[l,] <- esti_Z_gen[which(var_names == l),]
    } else {
      esti_Z_gen_mod[l,] <- pmin(1, colSums(esti_Z_gen[which(var_names == l),]))
    }
  }

  esti_Beta_gen <- matrix(nrow = q, ncol = p)
  diag(esti_Beta_gen) <- 0

  for(s in 1:p){
    if(categories[s] == 1){
      esti_Beta_gen[,s] <- esti_Beta[,which(var_names == s)]
    } else {
      esti_Beta_gen[,s] <- rowSums(esti_Beta[,which(var_names == s)])
    }
  }

  esti_Beta_gen_mod <- matrix(nrow = p, ncol = p)
  for(l in 1:p){
    if(categories[l ] == 1){
      esti_Beta_gen_mod[l,] <- esti_Beta_gen[which(var_names == l),]
    } else {
      esti_Beta_gen_mod[l,] <- colSums(esti_Beta_gen[which(var_names == l),])
    }
  }

  esti_Beta_gen_mod <- esti_Beta_gen_mod*lower.tri(esti_Beta_gen_mod)
  esti_Beta_gen_mod <- esti_Beta_gen_mod + t(esti_Beta_gen_mod)

  return(list("Adj_Beta" = esti_Beta_gen_mod, "Adj_Z" = esti_Z_gen_mod))
}
