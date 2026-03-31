# Implementation Plan: Enhancing Color Palettes and Removing `split.by`

This document details the plan to fix broken color palettes, add specialized palette categories, improve "Manual" color inputs, and remove the `split.by` menu configuration for Violin and Dot plots.

## User Review Required

> [!IMPORTANT]
> Please review the proposed changes below. Once approved, I will implement the changes and ensure this plan is kept version-controlled in the `docs` folder.

## Proposed Changes

### `app.R`
#### [MODIFY] `app.R` 
1. **Remove `split.by` Option for ViolinPlot and DotPlot** 
   - Update the dynamic UI generating script so the `split_by` dropdown evaluates whether `ptype %in% c("DimPlot", "FeaturePlot")`. `ViolinPlot` and `DotPlot` will no longer show it.
   - Remove `split_by` mapping dependencies on the server side for Violin plots (`DotPlot` already doesn't use it, but we will ensure it remains cleanly separated in UI).

2. **Fix Color Palettes & Support Categories**
   - *Current bug identified:* The UI exposes `palette_choice`, while the `color_utils.R` engine expects `palette_name` (returning `NULL` and resulting in default R colors silently). I will correct this misnomer.
   - Then, replace the palette dropdown options with categorized dropdowns exactly as requested: 
     - **Continuous:** `viridis`, `plasma`, `magma`, `inferno`, `cividis`, `Blues`, `Reds`, etc.
     - **Diverging:** `RdBu`, `RdYlBu`, `RdYlGn`, `Spectral`, etc.
     - **Discrete:** `Set1`, `Set2`, `Set3`, `Dark2`, `Paired`, etc.

3. **Improve the "Manual" Color Method**
   - Provide a sub-UI for "Manual Coloring" that lets users pick between:
     - **Color Pickers**: A visual pop-up window utilizing `colourpicker::colourInput` for each available level in the group.
     - **Hex Codes**: A text area where users can type/paste a comma-separated list of hex codes corresponding to the groups.

### `color_utils.R`
#### [MODIFY] `color_utils.R`
- Extend the `get_plot_colors` functionality to correctly parse comma-separated Hex codes.
- Map the new palette names from the categories directly, removing any pseudo-dependencies. Re-factor palette fallbacks.

### `plot_violin.R`
#### [MODIFY] `plot_violin.R`
- Remove the `split_by` function argument to fully disallow it from the backend violin plotting script.

## Open Questions

- Should we set a "fallback" continuous scheme (like Viridis) if someone enters fewer Hex codes than there are groups in the manual hex input mode? (My suggestion is to interpolate available colors or fallback to the standard gray if not enough are given).

## Verification Plan

### Manual Verification
1. Open the Docker container and load the app.
2. Verify `split.by` disappears appropriately on ViolinPlot and DotPlot selection.
3. Test that palettes dynamically update the plots on categorical data correctly.
4. Verify custom "Hex Code" inputs override values smoothly.
