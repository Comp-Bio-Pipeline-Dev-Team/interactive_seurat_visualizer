# ============================================================================
# Cluster Distribution Bar Plot Module
# ============================================================================
# This module provides a clean, self-contained function for generating
# cluster distribution bar plots (stacked bar charts showing proportions).
# ============================================================================

#' Create a cluster distribution bar plot
#'
#' @param obj Seurat object
#' @param group1 Name of metadata column for x-axis (samples/conditions)
#' @param group2 Name of metadata column for fill (clusters/cell types)
#' @param colors Optional vector of colors for fill aesthetic
#' @param flip Logical, whether to flip coordinates (horizontal bars)
#' @param show_counts Logical, whether to show total cell counts above bars
#' @param count_size Numeric, size of count labels (default: 3)
#' @param title Optional plot title
#' @param title_size Numeric, title font size (default: 16)
#' @param show_legend Logical, whether to show legend (default: TRUE)
#'
#' @return ggplot object
#' @export
plot_cluster_distribution <- function(obj, 
                                     group1, 
                                     group2, 
                                     colors = NULL,
                                     flip = FALSE,
                                     show_counts = FALSE,
                                     count_size = 3,
                                     title = NULL,
                                     title_size = 16,
                                     show_legend = TRUE) {
  
  # Validate inputs
  if (!group1 %in% names(obj@meta.data)) {
    stop(paste("group1 variable", group1, "not found in metadata"))
  }
  if (!group2 %in% names(obj@meta.data)) {
    stop(paste("group2 variable", group2, "not found in metadata"))
  }
  
  # Create frequency table
  df <- as.data.frame(table(obj@meta.data[[group1]], obj@meta.data[[group2]]))
  colnames(df) <- c("Sample", "Cluster", "Count")
  
  # Create base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Count, fill = Cluster)) + 
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(y = "Proportion", x = group1, fill = group2) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  
  # Add cell counts if requested
  if (show_counts) {
    df_totals <- aggregate(Count ~ Sample, df, sum)
    colnames(df_totals) <- c("Sample", "Total")
    df_totals$label_y <- 1.05  # Position above bar
    
    p <- p + ggplot2::geom_text(
      data = df_totals, 
      ggplot2::aes(x = Sample, y = label_y, label = paste0("n=", Total)), 
      inherit.aes = FALSE, 
      size = count_size, 
      fontface = "bold"
    )
  }
  
  # Flip coordinates if requested
  if (flip) {
    p <- p + ggplot2::coord_flip()
  }
  
  # Add title if provided
  if (!is.null(title) && nchar(title) > 0) {
    p <- p + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title_size, hjust = 0.5))
  }
  
  # Handle legend
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::theme(legend.position = "right")
  }
  
  return(p)
}


#' Validate cluster distribution plot parameters from Shiny inputs
#'
#' @param input Shiny input object
#' @param ns Namespace function
#' @param obj Seurat object
#'
#' @return List with validated parameters, or NULL if validation fails
#' @export
validate_cluster_distribution_params <- function(input, ns, obj) {
  group1 <- input[[ns("cdb_group1")]]
  group2 <- input[[ns("cdb_group2")]]
  
  # Check required inputs
  if (is.null(group1) || is.null(group2)) {
    return(NULL)
  }
  
  # Check that variables exist in metadata
  if (!group1 %in% names(obj@meta.data) || !group2 %in% names(obj@meta.data)) {
    return(NULL)
  }
  
  # Return validated parameters
  return(list(
    group1 = group1,
    group2 = group2,
    flip = input[[ns("flip_coords")]],
    show_counts = input[[ns("show_counts")]],
    count_size = input[[ns("count_size")]],
    title = input[[ns("custom_title")]],
    title_size = input[[ns("title_size")]],
    show_legend = input[[ns("show_legend")]]
  ))
}
