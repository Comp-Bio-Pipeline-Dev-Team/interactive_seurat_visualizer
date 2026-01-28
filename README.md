# Seurat Shiny App

An interactive R Shiny application for visualizing Seurat objects (v3, v4, v5).

## Features
- Upload Seurat (`.rds`) files (up to 10GB)
- 4-plot grid with independent settings
- Interactive visualization: DimPlot, FeaturePlot, ViolinPlot, DotPlot, ClusterDistrBar
- UCell signature scoring
- Advanced color controls (Viridis, RColorBrewer, MetBrewer palettes)
- Data subsetting and filtering
- Export plots as PNG/PDF/JPG

## Quick Start

### Option 1: Docker (Recommended for Deployment)

```bash
docker-compose up -d
```

Then open: http://localhost:3838

See [DOCKER_INSTRUCTIONS.md](DOCKER_INSTRUCTIONS.md) for details.

### Option 2: Local R Installation

```r
source("setup.R")  # Install dependencies
library(shiny)
runApp(".")
```

## Prerequisites
- R (>= 4.0.0) OR Docker
- For R installation: Seurat, Shiny, ggplot2, and other packages (see setup.R)

## Usage

1. Launch the app
2. Upload your `.rds` file containing the Seurat object
3. Configure plots using the sidebar tabs
4. Click on any plot to jump to its settings
