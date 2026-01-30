# ============================================================================
# Color Utilities for Seurat Interactive Visualizer
# ============================================================================
# This module provides unified interfaces for color palette generation,
# reducing dependency on multiple conflicting libraries and simplifying
# color management throughout the application.
# ============================================================================

#' Get colors from any palette library with a unified interface
#'
#' @param palette_name Name of the palette (e.g., "viridis", "Set1", "Hiroshige")
#' @param n_colors Number of colors needed
#' @param named Logical, whether to return a named vector
#' @param level_names Character vector of level names for naming colors
#'
#' @return Character vector of hex colors
#' @export
get_palette_colors <- function(palette_name, n_colors, named = FALSE, level_names = NULL) {
  if (n_colors == 0) return(NULL)
  
  colors <- NULL
  
  # Viridis family
  if (palette_name %in% c("viridis", "magma", "plasma", "inferno", "cividis")) {
    option_map <- list(
      viridis = "D",
      magma = "A", 
      plasma = "C",
      inferno = "B",
      cividis = "E"
    )
    colors <- viridis::viridis(n_colors, option = option_map[[palette_name]])
  }
  # MetBrewer palettes
  else if (requireNamespace("MetBrewer", quietly = TRUE) && 
           palette_name %in% names(MetBrewer::MetPalettes)) {
    colors <- MetBrewer::met.brewer(palette_name, n_colors)
  }
  # RColorBrewer palettes
  else if (palette_name %in% c("Set1", "Set2", "Set3", "Dark2", "Paired", 
                               "Pastel1", "Pastel2", "Accent", "Spectral",
                               "RdYlBu", "RdYlGn", "Blues", "Greens", "Reds")) {
    max_n <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
    if (n_colors <= max_n) {
      colors <- RColorBrewer::brewer.pal(n_colors, palette_name)
    } else {
      # Interpolate if we need more colors than available
      colors <- colorRampPalette(RColorBrewer::brewer.pal(max_n, palette_name))(n_colors)
    }
  }
  
  # If palette not found, return NULL
  if (is.null(colors)) {
    warning(paste("Palette", palette_name, "not found. Returning NULL."))
    return(NULL)
  }
  
  # Add names if requested
  if (named && !is.null(level_names) && length(level_names) == n_colors) {
    names(colors) <- level_names
  }
  
  return(colors)
}


#' Get available palette names grouped by library
#'
#' @return Named list of palette names by library
#' @export
get_available_palettes <- function() {
  palettes <- list(
    viridis = c("viridis", "magma", "plasma", "inferno", "cividis"),
    rcolorbrewer = c("Set1", "Set2", "Set3", "Dark2", "Paired", "Pastel1", "Accent"),
    metbrewer = NULL
  )
  
  # Add MetBrewer palettes if available
  if (requireNamespace("MetBrewer", quietly = TRUE)) {
    palettes$metbrewer <- names(MetBrewer::MetPalettes)
  }
  
  return(palettes)
}


#' Get all available palette names as a flat vector
#'
#' @return Character vector of all palette names
#' @export
get_all_palette_names <- function() {
  all_palettes <- get_available_palettes()
  return(unlist(all_palettes, use.names = FALSE))
}


#' Extract manual colors from Shiny inputs
#'
#' @param input Shiny input object
#' @param ns Namespace function
#' @param levels Character vector of level names
#' @param named Logical, whether to return a named vector
#'
#' @return Character vector of hex colors
#' @export
get_manual_colors_from_inputs <- function(input, ns, levels, named = FALSE) {
  if (length(levels) == 0) return(NULL)
  
  colors <- sapply(levels, function(level) {
    color_input_id <- ns(paste0("col_", level))
    color_value <- input[[color_input_id]]
    
    # Return the color if available, otherwise default to gray
    if (is.null(color_value)) "gray" else color_value
  })
  
  # Add names if requested
  if (named) {
    names(colors) <- levels
  } else {
    # Remove names if not requested
    colors <- unname(colors)
  }
  
  return(colors)
}


#' Get colors for a plot based on source (Default, Palette, or Manual)
#'
#' @param input Shiny input object
#' @param ns Namespace function
#' @param levels Character vector of level names to color
#' @param named Logical, whether to return a named vector (for SCpubr compatibility)
#'
#' @return Character vector of hex colors, or NULL for default
#' @export
get_plot_colors <- function(input, ns, levels, named = FALSE) {
  color_source <- input[[ns("color_source")]]
  
  # Default colors - let plotting library handle it
  if (is.null(color_source) || color_source == "Default") {
    return(NULL)
  }
  
  n_colors <- length(levels)
  if (n_colors == 0) return(NULL)
  
  # Palette colors
  if (color_source == "Palette") {
    palette_name <- input[[ns("palette_name")]]
    if (is.null(palette_name)) return(NULL)
    
    return(get_palette_colors(palette_name, n_colors, named = named, level_names = levels))
  }
  
  # Manual colors
  if (color_source == "Manual") {
    return(get_manual_colors_from_inputs(input, ns, levels, named = named))
  }
  
  return(NULL)
}


#' Determine the grouping variable for color assignment based on plot type
#'
#' @param input Shiny input object
#' @param ns Namespace function
#' @param plot_type Type of plot being generated
#'
#' @return Character string of the grouping variable name
#' @export
get_color_grouping_variable <- function(input, ns, plot_type) {
  # ClusterDistrBar uses cdb_group2 (fill variable)
  if (plot_type == "ClusterDistrBar") {
    return(input[[ns("cdb_group2")]])
  }
  
  # All other plots use group_by
  group_var <- input[[ns("group_by")]]
  
  # Handle "Default" or "None" values
  if (!is.null(group_var) && group_var %in% c("Default", "None")) {
    return(NULL)
  }
  
  return(group_var)
}


#' Get levels for a grouping variable from a Seurat object
#'
#' @param obj Seurat object
#' @param group_var Name of grouping variable (metadata column)
#'
#' @return Character vector of unique levels
#' @export
get_group_levels <- function(obj, group_var) {
  # If no group variable specified, use active identity
  if (is.null(group_var) || group_var == "Default") {
    return(levels(Idents(obj)))
  }
  
  # Get levels from metadata
  return(sort(unique(obj@meta.data[[group_var]])))
}


#' Apply colors to a ggplot object
#'
#' @param plot ggplot object
#' @param colors Character vector of colors
#' @param aesthetic Which aesthetic to apply colors to ("fill" or "color")
#'
#' @return Modified ggplot object
#' @export
apply_colors_to_plot <- function(plot, colors, aesthetic = "fill") {
  if (is.null(colors)) return(plot)
  
  if (aesthetic == "fill") {
    plot <- plot + ggplot2::scale_fill_manual(values = colors)
  } else if (aesthetic == "color") {
    plot <- plot + ggplot2::scale_color_manual(values = colors)
  }
  
  return(plot)
}
