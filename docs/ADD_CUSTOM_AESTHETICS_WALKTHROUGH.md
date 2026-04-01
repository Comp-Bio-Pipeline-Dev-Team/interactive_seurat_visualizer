# Walkthrough: Custom Aesthetics (Palette Reversal & Axis Titles)

As requested, we implemented two new globally available styling controls for all plots in the Seurat Interactive Visualizer to allow for precise aesthetic customizations: **Reverse Palette Ordering** and **Axis Label Sizing**.

## Summary of Changes

1. ✅ **Reverse Palette UI & Logic**: Added a dynamic checkbox `[ ] Reverse Palette Order` under the "Palette" dropdown. If clicked, the backend `color_utils.R` dynamically inverts the parsed hex-array output mathematically using `rev(colors)`. This perfectly inverts the coloring logic for *all submodules* automatically.
2. ✅ **Gradient Mapping Upgrade**: For Continuous data like FeaturePlots, the mapping originally just picked the "lowest" and "highest" color and ignored the inner palette stops (rendering visually unfaithful MetBrewer gradients). We upgraded `plot_feature.R` to natively process full multi-color gradient arrays using `ggplot2::scale_color_gradientn()` so reversing complex multi-dimensional palettes actually looks accurate and beautiful!
3. ✅ **Axis Title Scaling**: We implemented `Axis Title Size` (default: 14) as a numeric input alongside the pre-existing *Axis Text Size* (now renamed *Axis Text Tick Size*). This globally intercepts `p <- p + theme(...)` for every generated plot, identically replicating how `Title Size` is scaled!

## Files Updated
- **`app.R`**: Added the UI components & hooked the new `axis_title_sz` variable into the centralized `generate_plot` ggplot2 theme wrapper. Removed color array truncation.
- **`color_utils.R`**: Added `reverse` parameter validation and `rev()` mapping to `get_palette_colors()`.
- **`plot_feature.R`**: Defended standard `scale_color_gradient` via strict `scale_color_gradientn` to protect dimensional fidelity on reversal.

## Usage
No new dependencies were added. Rebuild your environment (`docker build -t seurat-shiny-app .`) and spin it up to play with the settings! They appear natively under the `Appearance & Colors` tab for every single plot module.
