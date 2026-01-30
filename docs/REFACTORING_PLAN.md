# Refactor App for Reduced Complexity

## Goal

Refactor the Rshiny application to reduce complexity by:
1. Creating unified interfaces for color palettes
2. Extracting plot-specific logic into modular functions
3. Reducing scattered conditional logic
4. Making it easier to maintain and debug

## Current Problems

1. **Color palette complexity**: Three different libraries (viridis, MetBrewer, RColorBrewer) with different APIs
2. **Plot generation complexity**: One massive `generate_plot()` function with nested conditionals
3. **Library conflicts**: SCpubr vs Seurat plotting with different parameter names
4. **Maintenance burden**: Changes require touching multiple scattered locations

## Proposed Architecture

### Phase 1: Unified Color System

Create `color_utils.R` with:

- **`get_palette_colors(palette_name, n_colors)`**: Unified wrapper for all palette libraries
- **`get_manual_colors_from_inputs(input, ns, levels)`**: Extract manual color logic
- **`apply_colors_to_plot(plot, colors, aesthetic = "fill")`**: Consistent color application

**Benefits:**
- Single function call regardless of palette source
- Easy to add new palette libraries
- Centralized color logic

### Phase 2: Plot-Specific Modules

Create separate files for each plot type:

#### `plot_cluster_distribution.R`
```r
plot_cluster_distribution(obj, group1, group2, colors = NULL, 
                         flip = FALSE, show_counts = FALSE, ...)
```

#### `plot_dimension_reduction.R`
```r
plot_dimension_reduction(obj, reduction, group_by = NULL, 
                        split_by = NULL, colors = NULL, 
                        style = "Standard", ...)
```

#### `plot_feature.R`
```r
plot_feature(obj, features, reduction, split_by = NULL, 
            colors = NULL, style = "Standard", ...)
```

**Benefits:**
- Each plot type is self-contained
- Easy to test individually
- Clear parameter requirements
- Reduced cognitive load

### Phase 3: Simplified Main App

Refactor `generate_plot()` to:
```r
generate_plot <- function(id, obj) {
  ptype <- input[[ns("plot_type")]]
  
  # Get colors once, correctly
  colors <- get_plot_colors(id, obj, ptype)
  
  # Dispatch to appropriate plot function
  switch(ptype,
    "ClusterDistrBar" = plot_cluster_distribution(...),
    "DimPlot" = plot_dimension_reduction(...),
    "FeaturePlot" = plot_feature(...),
    ...
  )
}
```

## Implementation Order

### Step 1: Create Color Utilities
- [NEW] [color_utils.R](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/color_utils.R)

### Step 2: Create ClusterDistrBar Module (Proof of Concept)
- [NEW] [plot_cluster_distribution.R](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/plot_cluster_distribution.R)

### Step 3: Refactor app.R to Use New Modules
- [MODIFY] [app.R](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/app.R)
  - Source new utility files
  - Replace color logic with unified functions
  - Replace ClusterDistrBar logic with module function

### Step 4: Extend to Other Plot Types
- Create remaining plot modules
- Update app.R to use all modules

## Verification Plan

### Automated Tests
- Test color utilities with various inputs
- Test each plot module independently

### Manual Verification
- Load test Seurat object
- Test each plot type with:
  - Default colors
  - Palette colors (viridis, MetBrewer, RColorBrewer)
  - Manual colors
- Verify no visual regressions
