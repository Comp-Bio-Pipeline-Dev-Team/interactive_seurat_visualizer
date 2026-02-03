# renv Package Management

## Current Status

The project uses `renv` for reproducible package management. However, the `renv.lock` file currently only contains the renv package itself.

## To Update renv.lock

To capture all currently installed packages and their versions, run the following in R:

```r
# Activate renv
renv::activate()

# Snapshot current package state
renv::snapshot()
```

This will update `renv.lock` with all packages currently used by the application, including:
- shiny, shinythemes
- Seurat, SeuratObject
- ggplot2, patchwork, cowplot
- plotly, colourpicker
- MetBrewer, viridis, RColorBrewer
- DT, ggrepel
- UCell (optional)
- SCpubr (optional)
- And all their dependencies

## Docker Integration

The Dockerfile is configured to use renv:
1. Copies `renv.lock`, `.Rprofile`, and `renv/` files
2. Runs `renv::restore()` to install exact package versions
3. Ensures reproducible builds across environments

## Note

Until `renv::snapshot()` is run, the Dockerfile will fall back to installing packages directly. Once the lock file is populated, all builds will use the exact versions specified in `renv.lock`.
