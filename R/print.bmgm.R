#' Print Method for BMGM Objects
#'
#' @param x A bmgm object returned by \code{\link{bmgm}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the bmgm object.
#' @export
print.bmgm <- function(x, ...) {
  type <- x$type
  p <- length(type)
  n <- nrow(x$X)

  type_counts <- c(
    continuous = sum(type == "c"),
    discrete = sum(type == "d"),
    "zero-inflated" = sum(type == "z"),
    categorical = sum(type == "m")
  )
  type_counts <- type_counts[type_counts > 0]
  type_str <- paste(paste(type_counts, names(type_counts)), collapse = ", ")

  cat("Bayesian Mixed Graphical Model\n")
  cat(p, "variables (", type_str, ")\n")
  cat(n, "observations |", x$nburn, "burn-in |", x$nsample, "samples\n\n")

  adj <- x$adj_G
  edges <- which(adj != 0 & upper.tri(adj), arr.ind = TRUE)
  n_edges <- nrow(edges)

  cat("Edges detected:", n_edges, "\n")

  if (n_edges > 0) {
    var_labels <- if (!is.null(colnames(x$X))) colnames(x$X) else paste0("X", 1:p)
    for (i in seq_len(n_edges)) {
      cat("  ", var_labels[edges[i, 1]], "--", var_labels[edges[i, 2]], "\n")
    }
  }

  cat("\nUse $adj_G for the full adjacency matrix.\n")
  cat("Use $adj_Beta for estimated edge weights.\n")
  cat("Use $inclusion_probs for posterior edge inclusion probabilities.\n")
  cat("Use $bfdr_cutoff for the BFDR threshold used.\n")
  if (!is.null(x$adj_Beta_ce)) {
    cat("Use $adj_Beta_ce for context-specific (category-level) edge weights.\n")
  }

  invisible(x)
}
