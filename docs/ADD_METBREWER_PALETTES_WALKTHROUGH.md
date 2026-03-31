# Walkthrough: Adding MetBrewer Palettes

As requested, we have completed the original implementation: adding the `MetBrewer` palettes as a dedicated, standalone category in the main plot palette selector and the heatmap configuration.

## Summary

1. ✅ Added a single line to `app.R` modifying `plotControlUI` to map `MetBrewer::MetPalettes` into a new `MetBrewer` dropdown list under `palette_name`.
2. ✅ Restored the `hm_anno_palette` to group palettes by their parent library (`Viridis`, `RColorBrewer`, `MetBrewer`).
3. ✅ Retained the `color_utils.R` backend logic, which correctly renders `MetBrewer` colors when selected.
4. ✅ Kept all existing Continuous, Diverging, and Discrete palettes entirely functional and separated by category.

## Code Changes

### `app.R`
Added `"MetBrewer" = names(MetBrewer::MetPalettes)` to the palette group definitions under "Appearance & Colors". 

```diff
  conditionalPanel(
    condition = sprintf("input['%s'] == 'Palette'", ns("color_source")),
    selectInput(ns("palette_name"), "Palette", 
               choices = list(
                 "Continuous" = c("viridis", "plasma", "magma", "inferno", "cividis", "Blues", "Reds", "Greens"),
                 "Diverging" = c("RdBu", "RdYlBu", "RdYlGn", "Spectral"),
-                "Discrete" = c("Set1", "Set2", "Set3", "Dark2", "Paired", "Pastel1", "Pastel2", "Accent")
+                "Discrete" = c("Set1", "Set2", "Set3", "Dark2", "Paired", "Pastel1", "Pastel2", "Accent"),
+                "MetBrewer" = names(MetBrewer::MetPalettes)
               ))
  ),
```

---

## Testing Protocol & Docker Build Command

Because you are using `renv.lock` and Docker, all changes run cleanly without needing external dependencies installed on your Mac. 

### Build and Run

To apply this updated `app.R` without the Mac keychain blocking `sudo docker build`, remember to securely build the container without `sudo`:

```bash
# 1. Rebuild the image from the folder since app.R changed
docker build -t seurat-shiny-app .

# 2. Run the newly built image on port 3838
docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app
```

Then, visit `http://localhost:3838`.
1. Once loaded, click **Plot 1**'s `Appearance & Colors`.
2. Change the coloring method to **Palette**.
3. You will now see `MetBrewer` listed as an established category below `Continuous`, `Diverging`, and `Discrete`.
