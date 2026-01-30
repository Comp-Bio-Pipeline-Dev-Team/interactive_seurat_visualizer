# Implementation Plan: Modularize Remaining Plot Types (Updated)

## Goal

Complete the refactoring by extracting DimPlot, FeaturePlot, ViolinPlot, and DotPlot into self-contained modules using **custom ggplot2 functions** instead of relying on Seurat/SCpubr packages.

## Key Requirements

> [!IMPORTANT]
> **Custom Plotting Functions**
> - Write all plotting logic using ggplot2 directly
> - Do NOT rely on Seurat's DimPlot, FeaturePlot, VlnPlot, or DotPlot functions
> - Do NOT rely on SCpubr functions
> - Extract data from Seurat objects and plot manually

> [!IMPORTANT]
> **Assay and Layer Selection**
> - Default to `assay = "RNA"` and `layer = "data"`
> - Add parameters to allow users to select which assay and layer to use
> - Support: `layer = c("counts", "data", "scale.data")`

> [!IMPORTANT]
> **DotPlot Scaling**
> - Add parameter for expression type: `scale_expression = c("average", "z-score")`
> - Update legend to reflect which scaling is used
> - Default to average expression from "data" layer

## Proposed Modules

### 1. DimPlot Module (`plot_dimension_reduction.R`)

**Function signature:**
```r
plot_dimension_reduction(obj, 
                        reduction, 
                        group_by = NULL, 
                        split_by = NULL,
                        assay = "RNA",
                        colors = NULL, 
                        pt_size = 1,
                        title = NULL,
                        show_legend = TRUE,
                        flip = FALSE)
```

**Implementation approach:**
- Extract embedding coordinates from `obj@reductions[[reduction]]@cell.embeddings`
- Extract grouping from `obj@meta.data[[group_by]]`
- Build plot using `ggplot() + geom_point()`
- Handle splitting with `facet_wrap()` if needed

**Complexity:** Medium

---

### 2. FeaturePlot Module (`plot_feature.R`)

**Function signature:**
```r
plot_feature(obj, 
            features, 
            reduction, 
            split_by = NULL,
            assay = "RNA",
            layer = "data",
            colors = NULL, 
            pt_size = 1,
            title = NULL,
            show_legend = TRUE)
```

**Implementation approach:**
- Extract embedding coordinates from reduction
- Extract expression values from `GetAssayData(obj, assay = assay, layer = layer)[features, ]`
- Build plot using `ggplot() + geom_point() + scale_color_gradient()`
- Handle multiple features with `facet_wrap()`

**Complexity:** Medium

---

### 3. ViolinPlot Module (`plot_violin.R`)

**Function signature:**
```r
plot_violin(obj, 
           features, 
           group_by = NULL, 
           split_by = NULL,
           assay = "RNA",
           layer = "data",
           colors = NULL, 
           pt_size = 0.1,
           title = NULL,
           show_legend = TRUE,
           flip = FALSE)
```

**Implementation approach:**
- Extract expression from `GetAssayData(obj, assay = assay, layer = layer)[features, ]`
- Extract grouping from metadata
- Build plot using `ggplot() + geom_violin() + geom_jitter()`
- Handle splitting with `facet_wrap()` or dodge position

**Complexity:** Medium

---

### 4. DotPlot Module (`plot_dot.R`)

**Function signature:**
```r
plot_dot(obj, 
        features, 
        group_by = NULL,
        assay = "RNA",
        layer = "data",
        scale_expression = c("average", "z-score"),
        colors = NULL,
        title = NULL,
        show_legend = TRUE,
        flip = FALSE)
```

**Implementation approach:**
- Extract expression from `GetAssayData(obj, assay = assay, layer = layer)[features, ]`
- Calculate average expression per group
- Calculate percent expressed per group
- Apply z-scaling if `scale_expression = "z-score"`
- Build plot using `ggplot() + geom_point(aes(size = pct, color = avg_exp))`
- Update legend title based on scaling method

**Complexity:** High

---

## Implementation Strategy

### Phase 1: Create Modules (One at a time)

1. **DimPlot** - Most commonly used, simpler data extraction
2. **FeaturePlot** - Similar to DimPlot but with continuous scale
3. **ViolinPlot** - Different geometry but straightforward
4. **DotPlot** - Most complex (aggregation + scaling)

### Phase 2: Refactor app.R

Simplified dispatch pattern:

```r
generate_plot <- function(id, obj) {
  ns <- function(x) paste0(id, "-", x)
  ptype <- input[[ns("plot_type")]]
  req(ptype)
  
  # Common parameters
  assay <- input[[ns("assay")]] %||% "RNA"
  layer <- input[[ns("layer")]] %||% "data"
  color_group_var <- get_color_grouping_variable(input, ns, ptype)
  colors <- get_colors(id, obj, color_group_var)
  
  # Dispatch to module
  p <- switch(ptype,
    "ClusterDistrBar" = plot_cluster_distribution(...),
    "DimPlot" = plot_dimension_reduction(...),
    "FeaturePlot" = plot_feature(..., assay = assay, layer = layer),
    "ViolinPlot" = plot_violin(..., assay = assay, layer = layer),
    "DotPlot" = plot_dot(..., assay = assay, layer = layer)
  )
  
  return(p)
}
```

### Phase 3: Add UI Controls

Add to dynamic UI:
```r
selectInput(ns("assay"), "Assay", choices = names(obj@assays), selected = "RNA")
selectInput(ns("layer"), "Layer", choices = c("counts", "data", "scale.data"), selected = "data")
```

For DotPlot specifically:
```r
radioButtons(ns("dot_scale"), "Expression Scale", 
            choices = c("Average" = "average", "Z-score" = "z-score"),
            selected = "average")
```

### Phase 4: Testing Checklist

For each module, verify:
- ✅ Default colors work
- ✅ Palette colors work  
- ✅ Manual colors work
- ✅ **group_by parameter works correctly**
- ✅ **split_by parameter works correctly**
- ✅ **Different assays work (if multiple exist)**
- ✅ **Different layers work (counts, data, scale.data)**
- ✅ Flip coordinates works (where applicable)
- ✅ Custom titles work
- ✅ Legend toggle works

Additional for DotPlot:
- ✅ Average expression mode
- ✅ Z-score expression mode
- ✅ Legend reflects correct scaling

---

## Benefits

- **No package dependencies**: Not tied to Seurat/SCpubr API changes
- **Full control**: Can customize every aspect of plots
- **Flexibility**: Easy to add new features or modify behavior
- **Reduced complexity**: ~280 lines → ~50 lines in app.R (82% reduction)
- **User control**: Assay and layer selection for all plots

## Files to Create

- `plot_dimension_reduction.R` - Custom DimPlot using ggplot2
- `plot_feature.R` - Custom FeaturePlot using ggplot2
- `plot_violin.R` - Custom ViolinPlot using ggplot2
- `plot_dot.R` - Custom DotPlot using ggplot2

## Files to Modify

- `app.R` - Add assay/layer UI controls, simplify generate_plot()
- `Dockerfile` - Add new module files
- `docs/` - Update documentation with each commit

## Estimated Impact

| Metric | Current | After | Improvement |
|--------|---------|-------|-------------|
| generate_plot() lines | ~280 | ~50 | **82% reduction** |
| Package coupling | High | Low | **Independent** |
| User control | Limited | Full | **Assay/Layer selection** |
| Maintainability | Complex | Simple | **High** |
