# ============================================================================
# Dimension Reduction Plot Module
# ============================================================================
# This module provides custom ggplot2-based dimension reduction plotting
# (e.g., UMAP, tSNE, PCA) without relying on Seurat's DimPlot function.
# Compatible with Seurat v3, v4, and v5.
# ============================================================================

#' Plot dimension reduction (UMAP, tSNE, PCA, etc.)
#'
#' Creates a scatter plot of cells in reduced dimension space using custom
#' ggplot2 code. Extracts embedding coordinates directly from Seurat object
#' and builds the plot from scratch.
#'
#' @param obj Seurat object (v3, v4, or v5)
#'   Data type: S4 object of class "Seurat"
#'   
#' @param reduction Character string specifying which dimensional reduction to use
#'   Data type: character (length 1)
#'   Examples: "umap", "tsne", "pca"
#'   Must exist in obj@reductions
#'   
#' @param group_by Character string specifying metadata column for coloring cells
#'   Data type: character (length 1) or NULL
#'   Default: NULL (uses active identity)
#'   Must exist in obj@meta.data if specified
#'   
#' @param split_by Character string specifying metadata column for faceting
#'   Data type: character (length 1) or NULL
#'   Default: NULL (no faceting)
#'   Must exist in obj@meta.data if specified
#'   
#' @param colors Character vector of colors to use for groups
#'   Data type: character vector of hex colors or NULL
#'   Default: NULL (uses ggplot2 default colors)
#'   Length should match number of unique groups
#'   
#' @param pt_size Numeric value for point size
#'   Data type: numeric (length 1)
#'   Default: 1
#'   Range: typically 0.1 to 5
#'   
#' @param title Character string for plot title
#'   Data type: character (length 1) or NULL
#'   Default: NULL (no title)
#'   
#' @param show_legend Logical indicating whether to show legend
#'   Data type: logical (length 1)
#'   Default: TRUE
#'   
#' @param flip Logical indicating whether to flip coordinates
#'   Data type: logical (length 1)
#'   Default: FALSE
#'   
#' @return ggplot object that can be further customized or printed
#'   Data type: ggplot
#'
#' @examples
#' # Basic UMAP plot
#' p <- plot_dimension_reduction(seurat_obj, reduction = "umap")
#' 
#' # Colored by cell type with custom colors
#' p <- plot_dimension_reduction(seurat_obj, reduction = "umap",
#'                               group_by = "cell_type",
#'                               colors = c("#FF0000", "#00FF00", "#0000FF"))
#' 
#' # Split by condition
#' p <- plot_dimension_reduction(seurat_obj, reduction = "umap",
#'                               group_by = "cell_type",
#'                               split_by = "condition")
#'
#' @export
plot_dimension_reduction <- function(obj,
                                    reduction,
                                    group_by = NULL,
                                    split_by = NULL,
                                    colors = NULL,
                                    pt_size = 1,
                                    title = NULL,
                                    show_legend = TRUE,
                                    flip = FALSE) {
  
  # Validate inputs
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object (v3, v4, or v5)")
  }
  
  if (!reduction %in% names(obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in object. Available:",
               paste(names(obj@reductions), collapse = ", ")))
  }
  
  # Extract embedding coordinates
  # Works for Seurat v3, v4, and v5
  embeddings <- obj@reductions[[reduction]]@cell.embeddings
  
  # Get dimension names (e.g., "UMAP_1", "UMAP_2")
  dim_names <- colnames(embeddings)
  if (length(dim_names) < 2) {
    stop("Reduction must have at least 2 dimensions")
  }
  
  # Create data frame for plotting
  plot_df <- data.frame(
    dim1 = embeddings[, 1],
    dim2 = embeddings[, 2],
    cell = rownames(embeddings)
  )
  
  # Add grouping variable
  if (is.null(group_by)) {
    # Use active identity
    plot_df$group <- Idents(obj)
  } else {
    # Use specified metadata column
    if (!group_by %in% colnames(obj@meta.data)) {
      stop(paste("group_by column", group_by, "not found in metadata"))
    }
    plot_df$group <- obj@meta.data[[group_by]]
  }
  
  # Add splitting variable if specified
  if (!is.null(split_by)) {
    if (!split_by %in% colnames(obj@meta.data)) {
      stop(paste("split_by column", split_by, "not found in metadata"))
    }
    plot_df$split <- obj@meta.data[[split_by]]
  }
  
  # Create base plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = dim1, y = dim2, color = group)) +
    ggplot2::geom_point(size = pt_size, alpha = 0.8) +
    ggplot2::labs(
      x = dim_names[1],
      y = dim_names[2],
      color = if (is.null(group_by)) "Identity" else group_by
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }
  
  # Add faceting if split_by is specified
  if (!is.null(split_by)) {
    p <- p + ggplot2::facet_wrap(~ split, ncol = 2)
  }
  
  # Add title if specified
  if (!is.null(title) && nchar(title) > 0) {
    p <- p + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
  }
  
  # Handle legend visibility
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  # Flip coordinates if requested
  if (flip) {
    p <- p + ggplot2::coord_flip()
  }
  
  return(p)
}


#' Validate dimension reduction plot parameters from Shiny inputs
#'
#' Extracts and validates parameters from Shiny input object for creating
#' a dimension reduction plot. Returns NULL if required parameters are missing.
#'
#' @param input Shiny input object
#'   Data type: reactivevalues (Shiny input)
#'   
#' @param ns Namespace function for accessing namespaced inputs
#'   Data type: function
#'   Created with shiny::NS(id)
#'   
#' @param obj Seurat object to validate against
#'   Data type: S4 object of class "Seurat"
#'
#' @return Named list of validated parameters, or NULL if validation fails
#'   Data type: list or NULL
#'   List contains: reduction, group_by, split_by, pt_size, title, show_legend, flip
#'
#' @export
validate_dimension_reduction_params <- function(input, ns, obj) {
  reduction <- input[[ns("reduction")]]
  
  # Check required inputs
  if (is.null(reduction)) {
    return(NULL)
  }
  
  # Check that reduction exists
  if (!reduction %in% names(obj@reductions)) {
    return(NULL)
  }
  
  # Get optional parameters
  group_by <- input[[ns("group_by")]]
  if (!is.null(group_by) && group_by %in% c("Default", "None")) {
    group_by <- NULL
  }
  
  split_by <- input[[ns("split_by")]]
  if (!is.null(split_by) && split_by %in% c("None", "")) {
    split_by <- NULL
  }
  
  # Return validated parameters
  return(list(
    reduction = reduction,
    group_by = group_by,
    split_by = split_by,
    pt_size = input[[ns("pt_size")]],
    title = input[[ns("custom_title")]],
    show_legend = input[[ns("show_legend")]],
    flip = input[[ns("flip_coords")]]
  ))
}
