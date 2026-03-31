# Add MetBrewer Palettes Implementation Plan

Add the `MetBrewer` palettes to the main plot palette selector, ensuring existing continuous, diverging, and discrete palettes remain fully functional.

## User Review Required

Please review this simple plan, and upon approval, I will execute the changes, update the git repository, and ensure all planning documents are synced to the `docs/` folder.

## Proposed Changes

### `app.R`
#### [MODIFY] [app.R](file:///Volumes/USB-HDD/seurat_app/interactive_seurat_visualizer/app.R)

1. **Update `plotControlUI` Palette Dropdown**
   - Add a new "MetBrewer" group under the `"palette_name"` `selectInput` choices that dynamically pulls the list of MetBrewer palettes using `names(MetBrewer::MetPalettes)`.
   - The dropdown list will change from:
     ```r
     choices = list(
       "Continuous" = c(...),
       "Diverging" = c(...),
       "Discrete" = c(...)
     )
     ```
     To:
     ```r
     choices = list(
       "Continuous" = c(...),
       "Diverging" = c(...),
       "Discrete" = c(...),
       "MetBrewer" = names(MetBrewer::MetPalettes)
     )
     ```

*Note: The backend (`color_utils.R`) securely handles the `MetBrewer` color retrieval logic already, reverting smoothly if called. `MetBrewer` is already loaded globally (`library(MetBrewer)`) and exists within `renv.lock`.*

### `docs/`
#### [NEW] [docs/ADD_METBREWER_PALETTES_PLAN.md](file:///Volumes/USB-HDD/seurat_app/interactive_seurat_visualizer/docs/ADD_METBREWER_PALETTES_PLAN.md)
- Save this exact plan to the `docs/` folder as requested.

## Verification Plan

### Automated Tests
- Build Docker image successfully from the existing `renv.lock` if needed (the lockfile shouldn't change).

### Manual Verification
1. Open the app after implementation.
2. Select any plot (e.g., `DimPlot` or `FeaturePlot`).
3. Switch "Color Source" to `Palette`.
4. Open the `Palette` dropdown and verify the new `MetBrewer` sub-category exists.
5. Select a MetBrewer palette (e.g., "Cassatt1") and confirm the plot updates successfully.
6. Verify existing palettes under Continuous/Diverging/Discrete still properly update the plots.
