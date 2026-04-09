#' Plot Method for BMGM Objects
#'
#' @description Visualizes the estimated graphical model using \pkg{qgraph}.
#' Two display modes are available: edge weights (default) and edge inclusion
#' probabilities.
#'
#' @param x A bmgm object returned by \code{\link{bmgm}}.
#' @param mode Character string specifying the plot type. \code{"weights"}
#'   (default) displays signed edge weights from \code{adj_Beta}, with edges
#'   involving categorical variables shown in grey. \code{"inclusion"} displays
#'   posterior edge inclusion probabilities, with the BFDR cutoff shown in the
#'   title.
#' @param layout Character string specifying the graph layout. Default is
#'   \code{"circle"}. See \code{\link[qgraph]{qgraph}} for options.
#' @param labels Optional character vector of node labels. Defaults to column
#'   names of \code{X} or \code{X1, X2, ...}.
#' @param ... Additional arguments passed to \code{\link[qgraph]{qgraph}}.
#'
#' @return Invisibly returns the qgraph object.
#'
#' @details
#' In \code{"weights"} mode, positive edges are shown in green and negative
#' edges in red, following qgraph conventions. Edges involving at least one
#' categorical (\code{"m"}) variable are shown in grey. This is because the
#' collapsed edge weight sums across category-level associations, which may
#' have opposite signs (e.g., category 2 positively associated while category 3
#' is negatively associated). A grey edge indicates a significant conditional
#' dependence exists, but its direction is ambiguous at the variable level.
#' To inspect the signed category-level effects, use the context-specific graph
#' returned in \code{fit$adj_Beta_ce}.
#'
#' In \code{"inclusion"} mode, edges are shown in blue with width proportional
#' to the inclusion probability. Only edges present in the selected graph
#' (i.e., those that passed the BFDR threshold) are displayed.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 200
#' x1 <- rnorm(n)
#' X <- cbind(x1, rnorm(n, 0.5*x1), rnorm(n), sample(1:3, n, replace = TRUE))
#' type <- c("c", "c", "c", "m")
#' fit <- bmgm(X, type, nburn = 1000, nsample = 1000)
#' plot(fit)                    # Edge weights (default)
#' plot(fit, mode = "inclusion") # Inclusion probabilities
#' }
plot.bmgm <- function(x, mode = c("weights", "inclusion"), layout = "circle",
                       labels = NULL, ...) {

  if (!requireNamespace("qgraph", quietly = TRUE)) {
    stop("Package 'qgraph' is required for plotting. Install it with: install.packages('qgraph')")
  }

  mode <- match.arg(mode)
  type <- x$type
  p <- length(type)

  if (is.null(labels)) {
    cn <- colnames(x$X)
    if (is.null(cn) || all(cn == "") || all(is.na(cn))) {
      labels <- paste0("X", 1:p)
    } else {
      # Replace any blank names with X1, X2, etc.
      blank <- cn == "" | is.na(cn)
      cn[blank] <- paste0("X", which(blank))
      labels <- cn
    }
  }

  # Check if any edges exist
  adj_G <- x$adj_G
  if (sum(adj_G) == 0) {
    message("No edges detected at the current BFDR threshold. Nothing to plot.")
    return(invisible(NULL))
  }

  # Node colors: white background for all
  node_colors <- rep("white", p)

  # Node shapes by variable type
  # circle = continuous, square = discrete, diamond = zero-inflated, triangle = categorical
  type_shapes <- c("c" = "circle", "d" = "square", "z" = "diamond", "m" = "triangle")
  node_shapes <- type_shapes[type]

  if (mode == "weights") {

    adj <- x$adj_Beta

    # Identify edges involving categorical variables
    cat_idx <- which(type == "m")
    is_cat_edge <- matrix(FALSE, p, p)
    if (length(cat_idx) > 0) {
      is_cat_edge[cat_idx, ] <- TRUE
      is_cat_edge[, cat_idx] <- TRUE
    }

    # Build edge color matrix: green for positive, red for negative, grey for categorical
    edge_colors <- matrix("#FFFFFF", p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        if (adj_G[i, j] != 0) {
          if (is_cat_edge[i, j]) {
            edge_colors[i, j] <- "#888888"
          } else if (adj[i, j] > 0) {
            edge_colors[i, j] <- "#2CA02C"
          } else {
            edge_colors[i, j] <- "#D62728"
          }
        }
      }
    }

    # Layout with legend panel
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1.5))

    graph <- qgraph::qgraph(
      abs(adj) * adj_G,
      layout = layout,
      labels = labels,
      edge.color = edge_colors,
      color = node_colors,
      shape = node_shapes,
      title = "BMGM \u2014 Edge Weights",
      ...
    )

    # Draw legend
    graphics::par(mar = c(2, 0, 2, 0.5))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))

    # --- Edge colors ---
    y_start <- 0.92
    graphics::text(0.5, y_start, "Edge sign", font = 2, cex = 0.9, adj = 0.5)

    graphics::segments(0.05, y_start - 0.07, 0.25, y_start - 0.07, col = "#2CA02C", lwd = 2.5)
    graphics::text(0.32, y_start - 0.07, "Positive", cex = 0.8, adj = 0)

    graphics::segments(0.05, y_start - 0.14, 0.25, y_start - 0.14, col = "#D62728", lwd = 2.5)
    graphics::text(0.32, y_start - 0.14, "Negative", cex = 0.8, adj = 0)

    if (any(type == "m")) {
      graphics::segments(0.05, y_start - 0.21, 0.25, y_start - 0.21, col = "#888888", lwd = 2.5)
      graphics::text(0.32, y_start - 0.21, "Categorical", cex = 0.8, adj = 0)
    }

    # --- Node shapes ---
    # Only show types present in the data
    y_shapes <- if (any(type == "m")) y_start - 0.35 else y_start - 0.28
    graphics::text(0.5, y_shapes, "Variable type", font = 2, cex = 0.9, adj = 0.5)

    shape_info <- list(
      c = list(label = "Continuous", pch = 21),
      d = list(label = "Discrete", pch = 22),
      z = list(label = "Zero-inflated", pch = 23),
      m = list(label = "Categorical", pch = 24)
    )

    y_pos <- y_shapes - 0.07
    for (tp in c("c", "d", "z", "m")) {
      if (any(type == tp)) {
        graphics::points(0.15, y_pos, pch = shape_info[[tp]]$pch, cex = 1.8,
                         bg = "white", col = "black", lwd = 1.5)
        graphics::text(0.32, y_pos, shape_info[[tp]]$label, cex = 0.8, adj = 0)
        y_pos <- y_pos - 0.07
      }
    }

  } else {
    # Inclusion probability mode
    incl <- x$inclusion_probs
    cutoff <- x$bfdr_cutoff

    if (is.null(incl)) {
      stop("Inclusion probabilities not available in this fit object.")
    }

    # Only show edges that were selected
    incl_display <- incl * adj_G

    # Color gradient: light blue (low prob) to dark blue (high prob)
    blue_ramp <- grDevices::colorRampPalette(c("#BDD7E7", "#08519C"))
    n_colors <- 100
    color_palette <- blue_ramp(n_colors)

    # Scale color from cutoff to 1 (all displayed edges are above cutoff)
    edge_colors_incl <- matrix("#FFFFFF", p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        if (adj_G[i, j] != 0) {
          # Map [cutoff, 1] to [1, n_colors]
          scaled <- (incl[i, j] - cutoff) / max(1 - cutoff, 1e-10)
          idx <- max(1, ceiling(scaled * n_colors))
          edge_colors_incl[i, j] <- color_palette[idx]
        }
      }
    }

    # Use layout with extra margin on right for legend (narrow bar)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(6, 1))

    # Build edge labels with probabilities
    edge_labels <- matrix("", p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        if (adj_G[i, j] != 0 & i < j) {
          edge_labels[i, j] <- round(incl[i, j], 2)
        }
      }
    }

    graph <- qgraph::qgraph(
      incl_display,
      layout = layout,
      labels = labels,
      edge.color = edge_colors_incl,
      edge.width = 1.5,
      edge.labels = edge_labels,
      edge.label.cex = 1.0,
      edge.label.bg = "white",
      edge.label.color = "black",
      color = node_colors,
      shape = node_shapes,
      title = paste0("BMGM \u2014 Inclusion Probabilities (BFDR cutoff: ", round(cutoff, 3), ")"),
      maximum = 1,
      ...
    )

    # Draw gradient legend (scaled from cutoff to 1)
    graphics::par(mar = c(5, 1, 5, 3))
    legend_y <- seq(cutoff, 1, length.out = n_colors)
    graphics::image(1, legend_y, matrix(legend_y, nrow = 1),
                    col = color_palette, axes = FALSE, xlab = "", ylab = "")
    tick_at <- pretty(c(cutoff, 1), n = 5)
    tick_at <- tick_at[tick_at >= cutoff & tick_at <= 1]
    graphics::axis(4, at = tick_at, labels = round(tick_at, 2), las = 1)
    graphics::mtext("Inclusion\nProbability", side = 3, line = 1, cex = 0.8)
  }

  invisible(graph)
}
