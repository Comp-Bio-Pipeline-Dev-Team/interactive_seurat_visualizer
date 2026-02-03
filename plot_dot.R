# ============================================================================
# Dot Plot Module
# ============================================================================
# This module provides custom ggplot2-based dot plotting for feature
# expression without relying on Seurat's DotPlot function.
# Compatible with Seurat v3, v4, and v5.
# ============================================================================

#' Plot feature expression as dot plot
#'
#' Creates a dot plot showing average expression and percent expressed for
#' features across groups. Uses custom ggplot2 code to extract and aggregate
#' expression data from specified assay and layer.
#'
#' @param obj Seurat object (v3, v4, or v5)
#'   Data type: S4 object of class "Seurat"
#'   
#' @param features Character vector of feature names to plot
#'   Data type: character vector
#'   Features must exist in the specified assay
#'   
#' @param group_by Character string specifying metadata column for grouping
#'   Data type: character (length 1) or NULL
#'   Default: NULL (uses active identity)
#'   Must exist in obj@meta.data if specified
#'   
#' @param assay Character string specifying which assay to use
#'   Data type: character (length 1)
#'   Default: "RNA"
#'   Must exist in obj@assays
#'   
#' @param layer Character string specifying which data layer to use
#'   Data type: character (length 1)
#'   Default: "data"
#'   Options: "counts", "data", "scale.data"
#'   
#' @param scale_expression Character string specifying how to scale expression
#'   Data type: character (length 1)
#'   Default: "average"
#'   Options: "average" (raw average), "z-score" (z-scaled across groups)
#'   
#' @param colors Character vector of colors for expression gradient
#'   Data type: character vector of hex colors or NULL
#'   Default: NULL (uses blue to red gradient)
#'   If provided, should be length 2 (low, high)
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
#' # Basic dot plot
#' p <- plot_dot(seurat_obj, features = c("CD3D", "CD8A", "CD4"))
#' 
#' # Grouped by cell type with z-score scaling
#' p <- plot_dot(seurat_obj, features = c("CD3D", "CD8A"),
#'               group_by = "cell_type",
#'               scale_expression = "z-score")
#' 
#' # Use counts layer
#' p <- plot_dot(seurat_obj, features = c("CD3D", "CD8A"),
#'               layer = "counts")
#'
#' @export
plot_dot <- function(obj,
                    features,
                    group_by = NULL,
                    assay = "RNA",
                    layer = "data",
                    scale_expression = "average",
                    colors = NULL,
                    title = NULL,
                    show_legend = TRUE,
                    flip = FALSE) {
  
  # Validate inputs
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object (v3, v4, or v5)")
  }
  
  if (!assay %in% names(obj@assays)) {
    stop(paste("Assay", assay, "not found in object. Available:",
               paste(names(obj@assays), collapse = ", ")))
  }
  
  if (!scale_expression %in% c("average", "z-score")) {
    stop("scale_expression must be either 'average' or 'z-score'")
  }
  
  # Extract expression data
  tryCatch({
    expr_data <- Seurat::GetAssayData(obj, assay = assay, layer = layer)
  }, error = function(e) {
    stop(paste("Could not extract data from assay", assay, "layer", layer, ":", e$message))
  })
  
  # Check that features exist
  missing_features <- features[!features %in% rownames(expr_data)]
  if (length(missing_features) > 0) {
    warning(paste("Features not found in assay:", paste(missing_features, collapse = ", ")))
    features <- features[features %in% rownames(expr_data)]
  }
  
  if (length(features) == 0) {
    stop("No valid features found in the specified assay")
  }
  
  # Get grouping variable
  if (is.null(group_by)) {
    groups <- Idents(obj)
  } else {
    if (!group_by %in% colnames(obj@meta.data)) {
      stop(paste("group_by column", group_by, "not found in metadata"))
    }
    groups <- obj@meta.data[[group_by]]
  }
  
  # Calculate average expression and percent expressed for each feature-group combination
  plot_data_list <- lapply(features, function(feat) {
    expr_vec <- as.numeric(expr_data[feat, ])
    
    # Calculate metrics for each group
    group_stats <- lapply(unique(groups), function(grp) {
      cells_in_group <- which(groups == grp)
      expr_in_group <- expr_vec[cells_in_group]
      
      data.frame(
        feature = feat,
        group = as.character(grp),
        avg_exp = mean(expr_in_group, na.rm = TRUE),
        pct_exp = sum(expr_in_group > 0, na.rm = TRUE) / length(expr_in_group) * 100,
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, group_stats)
  })
  
  plot_df <- do.call(rbind, plot_data_list)
  
  # Apply scaling if requested
  if (scale_expression == "z-score") {
    # Z-score normalize average expression across groups for each feature
    plot_df <- do.call(rbind, lapply(split(plot_df, plot_df$feature), function(feat_df) {
      feat_df$scaled_exp <- scale(feat_df$avg_exp)[, 1]
      return(feat_df)
    }))
    
    exp_column <- "scaled_exp"
    legend_title <- "Z-score\nExpression"
  } else {
    exp_column <- "avg_exp"
    legend_title <- "Average\nExpression"
  }
  
  # Create dot plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes_string(x = "group", y = "feature", 
                                                     size = "pct_exp", color = exp_column)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(name = "Percent\nExpressed", range = c(0, 10)) +
    ggplot2::labs(
      x = if (is.null(group_by)) "Identity" else group_by,
      y = "Features",
      color = legend_title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )
  
  # Apply color scale
  if (!is.null(colors) && length(colors) >= 2) {
    p <- p + ggplot2::scale_color_gradient(low = colors[1], high = colors[length(colors)],
                                           name = legend_title)
  } else {
    # Default: blue (low) to red (high)
    p <- p + ggplot2::scale_color_gradient(low = "blue", high = "red",
                                           name = legend_title)
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


#' Validate dot plot parameters from Shiny inputs
#'
#' Extracts and validates parameters from Shiny input object for creating
#' a dot plot. Returns NULL if required parameters are missing.
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
#'   List contains: features, group_by, assay, layer, scale_expression, title, show_legend, flip
#'
#' @export
validate_dot_plot_params <- function(input, ns, obj) {
  features <- input[[ns("feature")]]
  
  # Check required inputs
  if (is.null(features)) {
    return(NULL)
  }
  
  # Get assay and layer
  assay <- input[[ns("assay")]]
  if (is.null(assay)) assay <- "RNA"
  
  layer <- input[[ns("layer")]]
  if (is.null(layer)) layer <- "data"
  
  # Get scaling method
  scale_expression <- input[[ns("dot_scale")]]
  if (is.null(scale_expression)) scale_expression <- "average"
  
  # Get optional parameters
  group_by <- input[[ns("group_by")]]
  if (!is.null(group_by) && group_by %in% c("Default", "None")) {
    group_by <- NULL
  }
  
  # Return validated parameters
  return(list(
    features = features,
    group_by = group_by,
    assay = assay,
    layer = layer,
    scale_expression = scale_expression,
    title = input[[ns("custom_title")]],
    show_legend = input[[ns("show_legend")]],
    flip = input[[ns("flip_coords")]]
  ))
}
