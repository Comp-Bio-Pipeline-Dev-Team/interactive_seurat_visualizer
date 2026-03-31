# Seurat Interactive Visualizer

A comprehensive R Shiny application for interactive visualization and analysis of Seurat single-cell RNA-seq objects.

## Features

### Visualization
- **4-Plot Grid**: Independent settings for each plot with click-to-activate controls
- **Plot Types**: DimPlot, FeaturePlot, ViolinPlot, DotPlot, ClusterDistrBar
- **SCpubr Integration**: Publication-ready plots with customizable aesthetics
- **Color Palettes**: Viridis, RColorBrewer, MetBrewer, or manual color selection
- **Interactive Controls**: Point size, orientation, split.by, group.by

### File Support
- Seurat objects (`.rds`)

## Quick Start (Docker)

### Option 1: download pre-built docker container (recommended)

The app runs in Docker, with all R packages version controlled  via `renv.lock` for reproducibility.  In this option, there is no need to build the container.  User can used a pre-built container from dockerhub.

```bash
# 1. Download image from DockerHub
sudo docker pull tbrunetti/interactive_seuart_visualizer:dev_03302026

# 2. Run the container on the command line
sudo docker run -d -p 3838:3838 tbrunetti/interactive_seuart_visualizer:dev_03302026

# 3. Access the app by opening your favorite web browser and going to the following link
http://localhost:3838 **or** http://0.0.0.0:3838
```

### Option 2: build docker container and run app

The app runs in Docker, with all R packages pinned via `renv.lock` for reproducibility.

```bash
# 1. Build the image (restores all packages from renv.lock)
sudo docker build -t seurat-shiny-app .

# 2. Run the container
sudo docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app

# 3. Access the app
open http://localhost:3838
```

See [`docs/QUICKSTART.md`](docs/QUICKSTART.md) and [`docs/DOCKER_COMMANDS.md`](docs/DOCKER_COMMANDS.md) for full details.

## Usage Guide

### 1. Upload Your Data
- Click "Browse" in the sidebar
- Select your Seurat object (`.rds`, `.h5ad`, or `.h5seurat`)
- Wait for the upload and validation

### 2. Visualization Tab

#### Configure Plots
- Use the sidebar tabs: **Plot 1**, **Plot 2**, **Plot 3**, **Plot 4**
- Each plot has independent settings:
  - **Plot Type**: Choose from DimPlot, FeaturePlot, ViolinPlot, DotPlot, ClusterDistrBar
  - **Features**: Type to search genes or metadata (autocomplete enabled)
  - **Grouping**: Select metadata columns for grouping/splitting
  - **Colors**: Choose Default, Palette, or Manual colors
  - **Style**: Toggle between Standard Seurat and SCpubr styles

#### Interactive Features
- **Click any plot** to jump to its settings tab
- **Active plot** is highlighted with a bold border
- **Preview plot** in Subsetting tab shows your full dataset

### 3. Export Tab
- Choose target (single plot or all 4 plots)
- Select format (PNG, PDF, JPG)
- Set dimensions
- Download

## Prerequisites

- Docker (install from [docker.com](https://docker.com))
- All R packages are managed via `renv.lock` — no local R installation required

## Tips & Tricks

- **Large Datasets**: SCpubr automatically uses rasterization for >50,000 cells
- **Many Cell Types**: Manual color picker is scrollable for >20 levels
- **Feature Search**: Start typing in feature boxes for instant autocomplete
- **Split Plots**: Use `split.by` to facet plots by experimental conditions
- **Orientation**: Flip DotPlot/ViolinPlot to horizontal for better readability

## Troubleshooting

**App won't start**: Check that all dependencies are installed via `setup.R`

**Upload fails**: Ensure your file is a valid Seurat object and <10GB

**Plots not updating**: Click the plot to activate it and check its settings tab

**SCpubr errors**: Some metadata combinations may not be supported; try Standard style

## Project Structure

```
interactive_seurat_visualizer/
├── app.R                         # Main Shiny application
├── color_utils.R                 # Unified color palette system
├── plot_dimension_reduction.R    # DimPlot module
├── plot_feature.R                # FeaturePlot module
├── plot_violin.R                 # ViolinPlot module
├── plot_dot.R                    # DotPlot module
├── plot_cluster_distribution.R   # ClusterDistrBar module
├── ui_landing_page.R             # Landing page UI
├── server_landing_page.R         # Landing page server logic
├── enrichment_backend.R          # Gene enrichment analysis
├── setup.R                       # App initialization helpers
├── renv.lock                     # R package versions (source of truth for Docker)
├── Dockerfile                    # Docker image (builds from renv.lock)
├── docker-compose.yml            # Docker orchestration
└── docs/                         # All documentation (see docs/README.md)
```

See [`docs/README.md`](docs/README.md) for the full documentation index.

## Citation

If you use this app in your research, please cite:
- **Seurat**: Hao et al., Cell (2021)
- **SCpubr**: Blanco-Carmona, E. (2022)
- **UCell**: Andreatta & Carmona, Nat Comput Sci (2021)

## License

MIT License - see LICENSE file for details
