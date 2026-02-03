# ============================================================================
# Feature Plot Module
# ============================================================================
# This module provides custom ggplot2-based feature expression plotting
# in reduced dimension space without relying on Seurat's FeaturePlot function.
# Compatible with Seurat v3, v4, and v5.
# ============================================================================

#' Plot feature expression in dimension reduction space
#'
#' Creates a scatter plot showing expression of one or more features overlaid
#' on dimension reduction coordinates. Uses custom ggplot2 code to extract
#' expression data from specified assay and layer.
#'
#' @param obj Seurat object (v3, v4, or v5)
#'   Data type: S4 object of class "Seurat"
#'   
#' @param features Character vector of feature names to plot
#'   Data type: character vector
#'   Features must exist in the specified assay
#'   Multiple features will be faceted
#'   
#' @param reduction Character string specifying which dimensional reduction to use
#'   Data type: character (length 1)
#'   Examples: "umap", "tsne", "pca"
#'   Must exist in obj@reductions
#'   
#' @param split_by Character string specifying metadata column for faceting
#'   Data type: character (length 1) or NULL
#'   Default: NULL (no faceting by metadata)
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
#' @param colors Character vector of colors for gradient (low to high)
#'   Data type: character vector of hex colors or NULL
#'   Default: NULL (uses viridis color scale)
#'   If provided, should be length 2 (low, high)
#'   
#' @param pt_size Numeric value for point size
#'   Data type: numeric (length 1)
#'   Default: 1
#'   Range: typically 0.1 to 5
#'   
#' @param title Character string for plot title
#'   Data type: character (length 1) or NULL
#'   Default: NULL (uses feature name as title)
#'   
#' @param show_legend Logical indicating whether to show legend
#'   Data type: logical (length 1)
#'   Default: TRUE
#'
#' @return ggplot object that can be further customized or printed
#'   Data type: ggplot
#'
#' @examples
#' # Basic feature plot
#' p <- plot_feature(seurat_obj, features = "CD3D", reduction = "umap")
#' 
#' # Multiple features with custom colors
#' p <- plot_feature(seurat_obj, features = c("CD3D", "CD8A"),
#'                   reduction = "umap",
#'                   colors = c("lightgray", "red"))
#' 
#' # Use counts layer instead of data
#' p <- plot_feature(seurat_obj, features = "CD3D",
#'                   reduction = "umap", layer = "counts")
#'
#' @export
plot_feature <- function(obj,
                        features,
                        reduction,
                        split_by = NULL,
                        assay = "RNA",
                        layer = "data",
                        colors = NULL,
                        pt_size = 1,
                        title = NULL,
                        show_legend = TRUE) {
  
  # Validate inputs
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object (v3, v4, or v5)")
  }
  
  if (!assay %in% names(obj@assays)) {
    stop(paste("Assay", assay, "not found in object. Available:",
               paste(names(obj@assays), collapse = ", ")))
  }
  
  if (!reduction %in% names(obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in object. Available:",
               paste(names(obj@reductions), collapse = ", ")))
  }
  
  # Extract embedding coordinates
  embeddings <- obj@reductions[[reduction]]@cell.embeddings
  dim_names <- colnames(embeddings)
  
  if (length(dim_names) < 2) {
    stop("Reduction must have at least 2 dimensions")
  }
  
  # Extract expression data from specified assay and layer
  # Compatible with Seurat v3, v4, and v5
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
  
  # Create base data frame with coordinates
  plot_df <- data.frame(
    dim1 = embeddings[, 1],
    dim2 = embeddings[, 2],
    cell = rownames(embeddings)
  )
  
  # Add split variable if specified
  if (!is.null(split_by)) {
    if (!split_by %in% colnames(obj@meta.data)) {
      stop(paste("split_by column", split_by, "not found in metadata"))
    }
    plot_df$split <- obj@meta.data[[split_by]]
  }
  
  # Handle single vs multiple features
  if (length(features) == 1) {
    # Single feature - simple plot
    plot_df$expression <- as.numeric(expr_data[features, ])
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = dim1, y = dim2, color = expression)) +
      ggplot2::geom_point(size = pt_size, alpha = 0.8) +
      ggplot2::labs(
        x = dim_names[1],
        y = dim_names[2],
        color = "Expression"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.title = ggplot2::element_text(size = 11, face = "bold"),
        legend.text = ggplot2::element_text(size = 10)
      )
    
    # Apply color scale
    if (!is.null(colors) && length(colors) >= 2) {
      p <- p + ggplot2::scale_color_gradient(low = colors[1], high = colors[length(colors)])
    } else {
      p <- p + ggplot2::scale_color_viridis_c(option = "viridis")
    }
    
    # Add title
    if (is.null(title)) {
      title <- features
    }
    
  } else {
    # Multiple features - create long format and facet
    expr_list <- lapply(features, function(feat) {
      df <- plot_df
      df$expression <- as.numeric(expr_data[feat, ])
      df$feature <- feat
      return(df)
    })
    
    plot_df <- do.call(rbind, expr_list)
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = dim1, y = dim2, color = expression)) +
      ggplot2::geom_point(size = pt_size, alpha = 0.8) +
      ggplot2::facet_wrap(~ feature, ncol = 2) +
      ggplot2::labs(
        x = dim_names[1],
        y = dim_names[2],
        color = "Expression"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.title = ggplot2::element_text(size = 11, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        strip.text = ggplot2::element_text(size = 11, face = "bold")
      )
    
    # Apply color scale
    if (!is.null(colors) && length(colors) >= 2) {
      p <- p + ggplot2::scale_color_gradient(low = colors[1], high = colors[length(colors)])
    } else {
      p <- p + ggplot2::scale_color_viridis_c(option = "viridis")
    }
  }
  
  # Add split faceting if specified (in addition to feature faceting)
  if (!is.null(split_by)) {
    if (length(features) == 1) {
      p <- p + ggplot2::facet_wrap(~ split, ncol = 2)
    } else {
      p <- p + ggplot2::facet_grid(feature ~ split)
    }
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
  
  return(p)
}


#' Validate feature plot parameters from Shiny inputs
#'
#' Extracts and validates parameters from Shiny input object for creating
#' a feature plot. Returns NULL if required parameters are missing.
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
#'   List contains: features, reduction, split_by, assay, layer, pt_size, title, show_legend
#'
#' @export
validate_feature_plot_params <- function(input, ns, obj) {
  features <- input[[ns("feature")]]
  reduction <- input[[ns("reduction")]]
  
  # Check required inputs
  if (is.null(features) || is.null(reduction)) {
    return(NULL)
  }
  
  # Check that reduction exists
  if (!reduction %in% names(obj@reductions)) {
    return(NULL)
  }
  
  # Get assay and layer
  assay <- input[[ns("assay")]]
  if (is.null(assay)) assay <- "RNA"
  
  layer <- input[[ns("layer")]]
  if (is.null(layer)) layer <- "data"
  
  # Get optional parameters
  split_by <- input[[ns("split_by")]]
  if (!is.null(split_by) && split_by %in% c("None", "")) {
    split_by <- NULL
  }
  
  # Return validated parameters
  return(list(
    features = features,
    reduction = reduction,
    split_by = split_by,
    assay = assay,
    layer = layer,
    pt_size = input[[ns("pt_size")]],
    title = input[[ns("custom_title")]],
    show_legend = input[[ns("show_legend")]]
  ))
}
