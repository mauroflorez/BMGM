#' Plot Method for BMGM Objects
#'
#' @description Visualizes the estimated graphical model using \pkg{qgraph}.
#' Three display modes are available: edge weights (default), edge inclusion
#' probabilities, and context-specific graph.
#'
#' @param x A bmgm object returned by \code{\link{bmgm}}.
#' @param mode Character string specifying the plot type. \code{"weights"}
#'   (default) displays signed edge weights from \code{adj_Beta}, with edges
#'   involving categorical variables shown in grey. \code{"inclusion"} displays
#'   posterior edge inclusion probabilities, with the BFDR cutoff shown in the
#'   title. \code{"context"} displays the context-specific graph where each
#'   category of a multinomial variable is shown as a separate node, with
#'   signed edges (green/red) at the category level.
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
#' In \code{"context"} mode, the context-specific graph is displayed where
#' each category of a multinomial variable gets its own node (e.g., Group.2,
#' Group.3 for a 3-category variable, where category 1 is the reference).
#' This allows inspection of signed category-level effects that are hidden
#' by the grey edges in the default weights mode. Only available when the
#' model includes categorical variables.
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
plot.bmgm <- function(x, mode = c("weights", "inclusion", "context"), layout = "circle",
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

    # Hint about context-specific mode if categoricals are present
    if (any(type == "m")) {
      message("Note: Grey edges involve categorical variables (sign is ambiguous). ",
              "Use plot(fit, mode = \"context\") to see category-level signed effects.")
    }

  } else if (mode == "context") {
    # Context-specific mode: q x q graph with category-level nodes
    adj_ce <- x$adj_Beta_ce
    adj_Z_ce <- x$adj_Z_ce

    if (is.null(adj_ce) || is.null(adj_Z_ce)) {
      stop("Context-specific graph not available. Requires categorical variables and context_spec = TRUE.")
    }

    q <- ncol(adj_ce)

    # Use the tag names on adj_Beta_ce, or rebuild if missing
    ce_tags <- colnames(adj_ce)
    if (is.null(ce_tags)) {
      # Rebuild tags from type info
      categories <- sapply(1:p, function(s) {
        if (type[s] == "m") length(unique(x$X[!is.na(x$X[,s]), s])) else 1
      })
      ce_tags <- c()
      for (s in 1:p) {
        if (categories[s] == 1) {
          ce_tags <- c(ce_tags, paste0(s))
        } else {
          ce_tags <- c(ce_tags, paste0(s, "_", 1:categories[s]))
        }
      }
      if (length(ce_tags) != q) {
        stop("Cannot determine node labels for context-specific graph. Dimension mismatch.")
      }
    }

    cn <- colnames(x$X)
    if (is.null(cn) || all(cn == "") || all(is.na(cn))) {
      cn <- paste0("X", 1:p)
    } else {
      blank <- cn == "" | is.na(cn)
      cn[blank] <- paste0("X", which(blank))
    }

    # Map tags to readable labels: "3_2" -> "VarName.2"
    ce_labels <- sapply(ce_tags, function(tg) {
      parts <- strsplit(tg, "_")[[1]]
      var_idx <- as.integer(parts[1])
      if (length(parts) == 1) {
        cn[var_idx]
      } else {
        paste0(cn[var_idx], ".", parts[2])
      }
    }, USE.NAMES = FALSE)

    # Node shapes and colors
    ce_node_colors <- rep("white", q)
    ce_node_shapes <- character(q)
    for (k in seq_along(ce_tags)) {
      var_idx <- as.integer(strsplit(ce_tags[k], "_")[[1]][1])
      ce_node_shapes[k] <- unname(type_shapes[type[var_idx]])
    }

    # Build edge colors: green/red based on sign
    adj_ce_display <- adj_ce * (adj_Z_ce != 0)
    edge_colors_ce <- matrix("#FFFFFF", q, q)
    for (i in 1:q) {
      for (j in 1:q) {
        if (adj_Z_ce[i, j] != 0) {
          if (adj_ce[i, j] > 0) {
            edge_colors_ce[i, j] <- "#2CA02C"
          } else {
            edge_colors_ce[i, j] <- "#D62728"
          }
        }
      }
    }

    # Layout with legend
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1.5))

    graph <- qgraph::qgraph(
      abs(adj_ce_display),
      layout = layout,
      labels = ce_labels,
      edge.color = edge_colors_ce,
      color = ce_node_colors,
      shape = ce_node_shapes,
      title = "BMGM \u2014 Context-Specific Graph",
      ...
    )

    # Draw legend
    graphics::par(mar = c(2, 0, 2, 0.5))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))

    y_start <- 0.92
    graphics::text(0.5, y_start, "Edge sign", font = 2, cex = 0.9, adj = 0.5)

    graphics::segments(0.05, y_start - 0.07, 0.25, y_start - 0.07, col = "#2CA02C", lwd = 2.5)
    graphics::text(0.32, y_start - 0.07, "Positive", cex = 0.8, adj = 0)

    graphics::segments(0.05, y_start - 0.14, 0.25, y_start - 0.14, col = "#D62728", lwd = 2.5)
    graphics::text(0.32, y_start - 0.14, "Negative", cex = 0.8, adj = 0)

    y_shapes <- y_start - 0.28
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

    graphics::text(0.5, y_pos - 0.05, "Cat. nodes show\nindividual levels\n(e.g., Group.2)",
                   cex = 0.7, adj = 0.5, font = 3)

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
