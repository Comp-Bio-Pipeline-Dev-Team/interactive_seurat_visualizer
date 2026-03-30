# Seurat Interactive Visualizer

A comprehensive R Shiny application for interactive visualization and analysis of Seurat single-cell RNA-seq objects.

## Features

### Visualization
- **4-Plot Grid**: Independent settings for each plot with click-to-activate controls
- **Plot Types**: DimPlot, FeaturePlot, ViolinPlot, DotPlot, ClusterDistrBar
- **SCpubr Integration**: Publication-ready plots with customizable aesthetics
- **Color Palettes**: Viridis, RColorBrewer, MetBrewer, or manual color selection
- **Interactive Controls**: Point size, orientation, split.by, group.by

### Analysis
- **Differential Expression**: 
  - Group-vs-Group or One-vs-Rest comparisons
  - Interactive results table with sorting/filtering
  - Volcano plots with threshold lines
  - CSV export
- **UCell Signatures**: Calculate and visualize gene signature scores
- **Data Subsetting**: Filter cells by metadata or gene expression

### File Support
- Seurat objects (`.rds`)
- H5AD files (`.h5ad`)
- H5Seurat files (`.h5seurat`)

## Quick Start (Docker)

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

### 3. Differential Expression Tab

#### Run DE Analysis
1. Select **Comparison Group** (metadata column)
2. Choose **Ident 1** (group of interest)
3. Choose **Ident 2** (comparison group or "All Others")
4. Adjust thresholds (LogFC, Min Pct, Test Type)
5. Click **Run Differential Expression**

#### View Results
- **Results Table**: Sortable, searchable table with all markers
- **Volcano Plot**: Interactive visualization with threshold lines showing your cutoffs
- **Download CSV**: Export results for further analysis

### 4. Signatures Tab
- Enter genes (one per line)
- Name your signature
- Click **Calculate UCell**
- Signature scores appear in feature dropdowns as `[Name]_UCell`

### 5. Subsetting Tab
- **Filter by Metadata**: Select column and levels to keep
- **Preview**: See your data before filtering
- **Reset**: Restore original dataset

### 6. Export Tab
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
