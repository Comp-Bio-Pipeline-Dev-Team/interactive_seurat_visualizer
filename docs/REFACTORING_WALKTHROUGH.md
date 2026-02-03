# Refactoring Walkthrough: Modular Architecture

## Summary

Successfully refactored the Rshiny app to reduce complexity through modular design:
1. ✅ Fixed ClusterDistrBar manual color picker bug
2. ✅ Created unified color utility system
3. ✅ Extracted ClusterDistrBar into self-contained module
4. ✅ Reduced ~100 lines of scattered logic to ~30 lines of clean code

---

## Changes Made

### 1. Fixed ClusterDistrBar Color Picker Bug

**Issue:** Manual colors weren't applied because the wrong grouping variable was used.

**Fix in [app.R](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/app.R#L1332-L1344):**
```r
# Determine the correct grouping variable for color retrieval
# For ClusterDistrBar, use cdb_group2 (fill variable) instead of group_by
color_group_var <- grp
if (ptype == "ClusterDistrBar") {
  color_group_var <- input[[ns("cdb_group2")]]
}
```

### 2. Created Unified Color System

**New file:** [color_utils.R](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/color_utils.R)

**Key functions:**
- `get_palette_colors()` - Unified wrapper for viridis, MetBrewer, RColorBrewer
- `get_manual_colors_from_inputs()` - Extract manual colors from Shiny inputs
- `get_plot_colors()` - Single function to get colors regardless of source
- `get_group_levels()` - Consistent level extraction from Seurat objects

**Before (60+ lines):**
```r
# Scattered logic for each palette library
if (pal %in% c("viridis", "magma", ...)) {
  option_map <- list(...)
  cols <- viridis::viridis(n, option=option_map[[pal]])
} else if (pal %in% names(MetBrewer::MetPalettes)) {
  cols <- MetBrewer::met.brewer(pal, n)
} else if (pal %in% c("Set1", "Set2", ...)) {
  # More complex logic...
}
```

**After (10 lines):**
```r
# Clean wrapper using new utilities
get_colors <- function(id, obj, group_var, for_scpubr = FALSE) {
  ns <- function(x) paste0(id, "-", x)
  levels <- get_group_levels(obj, group_var)
  colors <- get_plot_colors(input, ns, levels, named = for_scpubr)
  return(colors)
}
```

### 3. Created ClusterDistrBar Module

**New file:** [plot_cluster_distribution.R](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/plot_cluster_distribution.R)

**Key functions:**
- `plot_cluster_distribution()` - Self-contained plotting function
- `validate_cluster_distribution_params()` - Parameter validation

**Before (36 lines of inline plotting):**
```r
if (ptype == "ClusterDistrBar") {
  g1 <- input[[ns("cdb_group1")]]
  g2 <- input[[ns("cdb_group2")]]
  req(g1, g2)
  df <- as.data.frame(table(...))
  # ... 30+ more lines of ggplot logic
}
```

**After (19 lines with clean function call):**
```r
if (ptype == "ClusterDistrBar") {
  params <- validate_cluster_distribution_params(input, ns, obj)
  if (!is.null(params)) {
    p <- plot_cluster_distribution(
      obj = obj,
      group1 = params$group1,
      group2 = params$group2,
      colors = colors,
      flip = isTRUE(params$flip),
      show_counts = isTRUE(params$show_counts),
      # ... clean parameter passing
    )
  }
}
```

### 4. Updated Docker Configuration

**Modified:** [Dockerfile](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/Dockerfile#L39-L45)

Added new utility files to container:
```dockerfile
COPY color_utils.R /app/
COPY plot_cluster_distribution.R /app/
```

---

## Code Metrics

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Color logic (lines) | 60+ | 10 | **83% reduction** |
| ClusterDistrBar logic (lines) | 36 | 19 | **47% reduction** |
| Total complexity | High (scattered) | Low (modular) | **Maintainable** |
| Library dependencies | Tightly coupled | Abstracted | **Flexible** |

---

## Testing Instructions (Docker-based)

### 1. Build the Docker Image

```bash
cd /home/tonya/Downloads/github/interactive_seurat_visualizer
docker-compose build
```

### 2. Run the Application

```bash
docker-compose up
```

Access at: http://localhost:3838

### 3. Test ClusterDistrBar with Manual Colors

1. **Upload a Seurat object** (or use test_seurat.rds if available)
2. **Select Plot 1 settings:**
   - Plot Type: `ClusterDistrBar`
   - X Axis (Sample): Choose a metadata column
   - Fill (Cluster): Choose another metadata column
3. **Set Color Source:** `Manual`
4. **Change colors** using the color pickers
5. **Verify:** Colors should update in the plot immediately

### 4. Test Other Plot Types (Regression Testing)

Test that other plot types still work with manual colors:
- ✅ DimPlot with manual colors
- ✅ FeaturePlot with palette colors
- ✅ ViolinPlot with manual colors
- ✅ DotPlot with palette colors

### 5. Stop the Container

```bash
docker-compose down
```

---

## Benefits of Refactoring

### ✅ Reduced Complexity
- Single source of truth for color logic
- Self-contained plot modules
- Easier to understand and debug

### ✅ Improved Maintainability
- Changes to color palettes only need updates in one place
- Plot logic is isolated and testable
- Clear separation of concerns

### ✅ Better Extensibility
- Easy to add new palette libraries
- Simple to create new plot types
- Consistent interface across all plots

### ✅ Reduced Library Conflicts
- Abstraction layer protects against API changes
- Can swap libraries without touching main code
- Easier to handle version conflicts

---

## Next Steps (Optional)

To complete the full refactoring:

1. **Create remaining plot modules:**
   - `plot_dimension_reduction.R` for DimPlot
   - `plot_feature.R` for FeaturePlot
   - `plot_violin.R` for ViolinPlot
   - `plot_dot.R` for DotPlot

2. **Further simplify generate_plot():**
   - Use switch statement for plot dispatch
   - Remove remaining inline plotting logic

3. **Add unit tests:**
   - Test color utilities independently
   - Test plot modules with mock data
