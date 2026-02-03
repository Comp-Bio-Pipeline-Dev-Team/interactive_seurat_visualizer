# Package Usage Audit

## Packages Currently Loaded

### Core Required (KEEP)
- ✅ **shiny** - UI framework (used extensively)
- ✅ **shinythemes** - UI themes (used in navbarPage)
- ✅ **Seurat** - Core analysis (used extensively)
- ✅ **ggplot2** - Base plotting (used in all plot modules)
- ✅ **patchwork** - Plot composition (used for combining plots)
- ✅ **cowplot** - Plot themes (used for plot_grid, themes)
- ✅ **DT** - Data tables (used for metadata table, DE results)

### Color & Visualization (KEEP)
- ✅ **MetBrewer** - Color palettes (used in color_utils.R, heatmaps)
- ✅ **viridis** - Color palettes (used in color_utils.R, feature plots)
- ✅ **RColorBrewer** - Color palettes (used in color_utils.R, heatmaps)
- ✅ **plotly** - Interactive plots (used for volcano plot: ggplotly)
- ✅ **colourpicker** - Color picker UI (used for volcano plot colors, heatmap colors)
- ❌ **ggrepel** - Label positioning (library loaded but NOT USED in actual code)

### Optional/Conditional (REVIEW)
- ❌ **SCpubr** - Alternative plotting (USER WANTS TO REMOVE)
- ⚠️ **UCell** - Gene signature scoring (conditionally loaded, check if used)
- ⚠️ **clusterProfiler, enrichplot, fgsea, msigdbr, DOSE, igraph** - Enrichment analysis (conditionally loaded)

## Recommendations

### Remove Immediately
1. **SCpubr** - User requested removal
2. **ggrepel** - Loaded but never used (no `ggrepel::` or `geom_text_repel` calls)

### Keep (Actually Used)
- All color palette packages (MetBrewer, viridis, RColorBrewer) - actively used
- colourpicker - used for volcano plot and heatmap color selection
- plotly - used for interactive volcano plot
- All core packages

### Investigate Further
- UCell - need to check if gene signature scoring is actually used
- Enrichment packages - need to check if enrichment analysis tab is functional

## Actions Needed

1. Remove from `app.R`:
   - Line 13: `library(ggrepel)`
   - Lines 16-21: SCpubr conditional loading
   - Lines 82-91: SCpubr UI controls
   - Lines 1530-1538: SCpubr plot style options

2. Remove from `setup.R`:
   - ggrepel
   - SCpubr installation line

3. Remove from `Dockerfile`:
   - Line with SCpubr GitHub installation
   - Can potentially remove `remotes` package if only used for SCpubr

4. Check if `remotes` is used elsewhere before removing
