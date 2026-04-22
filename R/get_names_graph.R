#' Get Names for Graph
#'
#' @param p Number of original variables
#' @param q Number of transformed variables
#' @param categories Category sizes
#' @param var_names Variable names
#' @return A list of graph names and edge indicators
#' @export
get_names_graph <- function(p, q, categories, var_names) {
  tag <- NULL
  for (i in 1:p) {
    tag <- c(tag, if (categories[i] == 1) paste0(i) else paste0(i, "_", 1:categories[i]))
  }

  tag_Beta <- outer(tag, tag, paste, sep = "-")
  ind_noedge <- outer(var_names, var_names, "==")
  diag(ind_noedge) <- TRUE

  return(list("tags" = tag, "tag_Beta" = tag_Beta, "indicators_noedge" = !ind_noedge))
}
