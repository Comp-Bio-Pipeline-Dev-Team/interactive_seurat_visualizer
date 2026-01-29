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
- Seurat objects (`.rds`) up to 10GB
- H5AD files (`.h5ad`) with automatic conversion
- H5Seurat files (`.h5seurat`)

## Quick Start

### Option 1: Local R Installation

1. **Install Dependencies**:
```r
source("setup.R")
```

2. **Launch the App**:
```r
library(shiny)
runApp(".")
```

3. **Access**: The app will open automatically in your default browser at `http://127.0.0.1:XXXX`

### Option 2: Docker (Recommended for Deployment)

1. **Build and Run**:
```bash
docker-compose up -d
```

2. **Access**: Open `http://localhost:3838` in your browser

3. **Stop**:
```bash
docker-compose down
```

See [DOCKER_INSTRUCTIONS.md](DOCKER_INSTRUCTIONS.md) for detailed Docker deployment instructions.

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

### For Local Installation
- R >= 4.0.0
- Required packages (installed via `setup.R`):
  - Core: `shiny`, `Seurat`, `ggplot2`, `patchwork`
  - Visualization: `SCpubr`, `MetBrewer`, `viridis`, `RColorBrewer`
  - Analysis: `UCell`, `DT`, `ggrepel`
  - Utilities: `colourpicker`, `cowplot`, `plotly`

### For Docker
- Docker and Docker Compose

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
seurat_shiny_app/
├── app.R                    # Main application
├── setup.R                  # Dependency installer
├── Dockerfile               # Docker image definition
├── docker-compose.yml       # Docker orchestration
├── DOCKER_INSTRUCTIONS.md   # Docker deployment guide
└── README.md               # This file
```

## Citation

If you use this app in your research, please cite:
- **Seurat**: Hao et al., Cell (2021)
- **SCpubr**: Blanco-Carmona, E. (2022)
- **UCell**: Andreatta & Carmona, Nat Comput Sci (2021)

## License

MIT License - see LICENSE file for details
