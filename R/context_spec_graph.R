#' Compute Context-Specific Graph
#'
#' @param q Number of variables
#' @param post_Bta Posterior Beta samples
#' @param post_Z Posterior adjacency samples
#' @param tags Variable names
#' @param bfdr Bayesian False Discovery Rate
#' @return Estimated graph
#' @export
context_spec_graph <- function(q, post_Bta, post_Z, tags, bfdr) {
  esti_Beta <- matrix(rep(0, q*q), nrow = q, ncol = q)
  esti_Beta[upper.tri(esti_Beta)] <- colMeans(post_Bta)
  esti_Beta <- esti_Beta + t(esti_Beta)
  colnames(esti_Beta) <- tags
  rownames(esti_Beta) <- tags

  esti_Z <- matrix(rep(0, q*q), nrow = q, ncol = q)

  #FDR:
  post_inclusion <- colMeans(post_Z)
  fdr_c <- function(c){
    n_edges <- sum(post_inclusion > c)
    if(n_edges == 0) return(1.0)
    E_fdr <- sum((1 - post_inclusion)*(post_inclusion > c)) / n_edges
    return(E_fdr)
  }

  pos_c <- seq(0, 1, by = 0.01)
  expected_fdr <- sapply(pos_c, fdr_c)
  pos <- pos_c[expected_fdr < bfdr]
  cutoff <- if(length(pos) > 0) min(pos) else 1.0

  esti_Z[upper.tri(esti_Z)] <- colMeans(post_Z)
  esti_Z <- esti_Z + t(esti_Z)

  esti_Z <- esti_Z*(esti_Z > cutoff)
  colnames(esti_Z) <- tags
  rownames(esti_Z) <- tags

  esti_Beta <- esti_Beta*esti_Z

  # Build full inclusion probability matrix (q x q)
  incl_prob <- matrix(0, nrow = q, ncol = q)
  incl_prob[upper.tri(incl_prob)] <- post_inclusion
  incl_prob <- incl_prob + t(incl_prob)
  colnames(incl_prob) <- tags
  rownames(incl_prob) <- tags

  return(list(ce_esti_Beta = esti_Beta, ce_esti_Z = esti_Z,
              inclusion_probs = incl_prob, bfdr_cutoff = cutoff))
}
