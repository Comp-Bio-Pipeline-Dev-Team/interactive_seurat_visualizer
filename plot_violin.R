# ============================================================================
# Violin Plot Module
# ============================================================================
# This module provides custom ggplot2-based violin plotting for feature
# expression without relying on Seurat's VlnPlot function.
# Compatible with Seurat v3, v4, and v5.
# ============================================================================

#' Plot feature expression as violin plots
#'
#' Creates violin plots showing distribution of feature expression across
#' groups. Uses custom ggplot2 code to extract expression data from
#' specified assay and layer.
#'
#' @param obj Seurat object (v3, v4, or v5)
#'   Data type: S4 object of class "Seurat"
#'   
#' @param features Character vector of feature names to plot
#'   Data type: character vector
#'   Features must exist in the specified assay
#'   Multiple features will be faceted
#'   
#' @param group_by Character string specifying metadata column for grouping
#'   Data type: character (length 1) or NULL
#'   Default: NULL (uses active identity)
#'   Must exist in obj@meta.data if specified
#'   
#' @param split_by Character string specifying metadata column for splitting violins
#'   Data type: character (length 1) or NULL
#'   Default: NULL (no splitting)
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
#' @param colors Character vector of colors for groups
#'   Data type: character vector of hex colors or NULL
#'   Default: NULL (uses ggplot2 default colors)
#'   Length should match number of unique groups
#'   
#' @param pt_size Numeric value for point size in jitter
#'   Data type: numeric (length 1)
#'   Default: 0.1
#'   Set to 0 to hide points
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
#' # Basic violin plot
#' p <- plot_violin(seurat_obj, features = "CD3D")
#' 
#' # Grouped by cell type with custom colors
#' p <- plot_violin(seurat_obj, features = "CD3D",
#'                  group_by = "cell_type",
#'                  colors = c("#FF0000", "#00FF00"))
#' 
#' # Multiple features
#' p <- plot_violin(seurat_obj, features = c("CD3D", "CD8A"),
#'                  group_by = "cell_type")
#'
#' @export
plot_violin <- function(obj,
                       features,
                       group_by = NULL,
                       split_by = NULL,
                       assay = "RNA",
                       layer = "data",
                       colors = NULL,
                       pt_size = 0.1,
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
  
  # Create data frame for plotting
  plot_list <- lapply(features, function(feat) {
    df <- data.frame(
      expression = as.numeric(expr_data[feat, ]),
      group = groups,
      feature = feat
    )
    
    # Add split variable if specified
    if (!is.null(split_by)) {
      if (!split_by %in% colnames(obj@meta.data)) {
        stop(paste("split_by column", split_by, "not found in metadata"))
      }
      df$split <- obj@meta.data[[split_by]]
    }
    
    return(df)
  })
  
  plot_df <- do.call(rbind, plot_list)
  
  # Create base plot
  if (is.null(split_by)) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = expression, fill = group)) +
      ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.7)
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = expression, fill = split)) +
      ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.7, position = ggplot2::position_dodge(0.9))
  }
  
  # Add points if requested
  if (pt_size > 0) {
    if (is.null(split_by)) {
      p <- p + ggplot2::geom_jitter(height = 0, width = 0.2, size = pt_size, alpha = 0.3)
    } else {
      p <- p + ggplot2::geom_jitter(height = 0, width = 0.2, size = pt_size, alpha = 0.3,
                                    position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9))
    }
  }
  
  # Add theme and labels
  p <- p +
    ggplot2::labs(
      x = if (is.null(group_by)) "Identity" else group_by,
      y = "Expression",
      fill = if (is.null(split_by)) (if (is.null(group_by)) "Identity" else group_by) else split_by
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  
  # Add faceting for multiple features
  if (length(features) > 1) {
    p <- p + ggplot2::facet_wrap(~ feature, scales = "free_y", ncol = 2) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 11, face = "bold"))
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


#' Validate violin plot parameters from Shiny inputs
#'
#' Extracts and validates parameters from Shiny input object for creating
#' a violin plot. Returns NULL if required parameters are missing.
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
#'   List contains: features, group_by, split_by, assay, layer, pt_size, title, show_legend, flip
#'
#' @export
validate_violin_plot_params <- function(input, ns, obj) {
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
    features = features,
    group_by = group_by,
    split_by = split_by,
    assay = assay,
    layer = layer,
    pt_size = input[[ns("pt_size")]],
    title = input[[ns("custom_title")]],
    show_legend = input[[ns("show_legend")]],
    flip = input[[ns("flip_coords")]]
  ))
}
