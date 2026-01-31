# Seurat Interactive Visualizer

A comprehensive R Shiny application for interactive visualization and analysis of Seurat single-cell RNA-seq objects.

## Features

### Visualization
- **4-Plot Grid**: Independent settings for each plot with click-to-activate controls
- **Plot Types**: DimPlot, FeaturePlot, ViolinPlot, DotPlot, ClusterDistrBar
- **SCpubr Integration**: Publication-ready plots with customizable aesthetics
- **Color Palettes**: Viridis, RColorBrewer, MetBrewer, or manual color selection
- **Interactive Controls**: Point size, orientation, split.by, group.by

### Gene Programs
- **Fast NMF Analysis**: Optimized gene program discovery using `snmf/r` (Sparse NMF with Alternating Least Squares)
- **Factor Heatmaps**: Visualize top genes per factor in clustered heatmaps
- **Spatial Visualization**: Project factor scores onto UMAP/tSNE
- **Gene Rankings**: Export top contributing genes for each program

### Analysis
- **Differential Expression**: 
  - Group-vs-Group or One-vs-Rest comparisons
  - Interactive results table with sorting/filtering
  - Volcano plots with customizable thresholds
  - CSV export
- **Pathway Enrichment**:
  - ORA (Over-Representation Analysis) and GSEA (Gene Set Enrichment Analysis)
  - Multiple databases: MSigDB Hallmarks, GO, KEGG, Reactome
  - Gene direction filtering for ORA (up/down/all regulated)
  - Visualization: Dot plots, bar plots, network plots (ORA), enrichment curves (GSEA)
  - Custom gene lists or direct integration with DE results
- **UCell Signatures**: Calculate and visualize gene signature scores
- **Data Subsetting**: Filter cells by metadata or gene expression

### Heatmap
- **Expression Heatmaps**: Generate SCpubr-style heatmaps with multiple genes
- **Z-score Scaling**: Row-wise normalization or raw expression
- **Grouping**: Organize cells by any metadata column
- **Customization**: Color palettes, dimensions, font sizes, rasterization for large datasets

### Data Management
- **File Support**: Seurat objects (.rds), H5AD files (.h5ad), H5Seurat files (.h5seurat) up to 10GB
- **Ensembl Conversion**: Convert Ensembl Gene IDs to Gene Symbols
- **QC Metrics**: Interactive quality control with custom mitochondrial column selection and filtering sliders

### Multimodal (CITE-seq) & Gating
- **Protein Visualization**: FeaturePlots for ADT data with full aesthetic control
- **Co-expression**: Scatterplots of Protein vs Gene expression with density contours
- **Gating**: Interactive flow-cytometry style gating (lasso/box select) to subset cells based on co-expression patterns

### VDJ Repertoire Analysis
- **Clonotype Visualization**: Analyze clonal expansion and distribution
- **Diversity Metrics**: Calculate Shannon, Simpson, and InvSimpson indices
- **Integration**: Link TCR/BCR data with scRNA-seq clusters
- **Export**: Visualizations and data tables for repertoire metrics


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

### 3. Gene Programs Tab

#### Run NMF Analysis
1. Set **Number of Factors (k)** (how many gene programs to discover)
2. Click **Run NMF** (may take 1-2 minutes for large datasets)
3. View results in three tabs:
   - **Factor Heatmap**: Top genes for each factor shown as clustered heatmap
   - **Factor Scores**: UMAP/tSNE colored by selected factor's expression
   - **Top Genes**: Interactive table of gene rankings per factor

#### Tips
- Start with k=5-10 for exploratory analysis
- Higher k values find more granular programs but increase computation time
- Factor scores automatically appear in feature dropdowns for visualization

### 4. Differential Expression Tab

#### Run DE Analysis
1. Select **Comparison Group** (metadata column)
2. Choose **Ident 1** (group of interest)
3. Choose **Ident 2** (comparison group or "All Others")
4. Adjust thresholds (LogFC, Min Pct, Test Type)
5. Click **Run Differential Expression**

#### View Results
- **Results Table**: Sortable, searchable table with all markers
- **Volcano Plot**: Interactive visualization with customizable threshold lines
- **Download CSV**: Export results for further analysis

#### Customize Volcano Plot
- Adjust point size, transparency, and colors
- Modify threshold lines to match your cutoffs
- Export as PNG, PDF, or JPG

### 5. Pathway Enrichment Tab

#### Input Sources
- **From DE Results**: Uses genes from your differential expression analysis
- **Custom Gene List**: Paste your own gene list (one per line)

#### Species Selection
- **Human/Mouse Support**: Explicitly choose the target organism for accurate database mapping (`org.Hs.eg.db` / `org.Mm.eg.db`)

#### Choose Analysis Type
- **ORA (Over-Representation Analysis)**: Tests if specific pathways are over-represented in your gene set
  - Use gene direction filter to analyze upregulated, downregulated, or all significant genes
  - Visualizations: Dot plot, bar plot, network plot
- **GSEA (Gene Set Enrichment Analysis)**: Tests for coordinated up/down-regulation of pathways
  - Uses full ranked gene list (requires DE results)
  - Visualizations: Enrichment curves showing running enrichment scores

#### Select Database
- **MSigDB Hallmarks**: Curated gene sets representing well-defined biological processes
- **GO Biological Process**: Gene Ontology biological processes
- **GO Molecular Function**: Gene Ontology molecular functions
- **GO Cellular Component**: Gene Ontology cellular components
- **KEGG**: Kyoto Encyclopedia of Genes and Genomes pathways
- **Reactome**: Reactome pathway database

#### Run Analysis
1. Configure input source and analysis type
2. Select organism (Human or Mouse)
3. Choose database
4. Adjust parameters (p-value, q-value, gene set sizes)
5. Click **Run Enrichment Analysis**

#### Export Results
- Download results table as CSV
- Download visualizations (dot plot, bar plot, network plot, or GSEA curves)
- Customize plot dimensions and format (PNG, PDF, JPG)

### 6. Heatmap Tab

#### Generate Expression Heatmap
1. Select **Features** (genes or metadata to visualize)
   - Use autocomplete to search for genes
   - Select multiple features
2. Choose **Group Cells By** (metadata column for organizing columns)
3. Select **Scaling** option:
   - **None**: Show raw expression values
   - **Row Z-score**: Normalize each gene across cells
4. Customize appearance:
   - **Color Palette**: Viridis, Magma, Plasma, Spectral, etc.
   - **Annotation Colors**: Choose palette or manual colors for grouping variable
   - **Dimensions**: Set width and height
   - **Font Sizes**: Adjust gene and legend font sizes
   - **Gap Width**: Control spacing between groups
5. Click **Generate Heatmap**
6. Download as PNG

#### Tips
- Enable **Rasterize** for large datasets to reduce file size
- Row Z-score is recommended for comparing genes with different expression ranges
- SCpubr styling ensures publication-ready aesthetics

### 7. Signatures Tab
- Enter genes (one per line)
- Name your signature
- Click **Calculate UCell**
- Signature scores appear in feature dropdowns as `[Name]_UCell`

### 8. Subsetting Tab
- **Filter by Metadata**: Select column and levels to keep
- **Preview**: See your data before filtering
- **Reset**: Restore original dataset

### 9. Multimodal (ADT) Tab
- **Protein Expression**: Visualize surface protein levels
- **Co-expression**: Plot Protein vs Gene expression
- **Aesthetics**: Full control over point size, palettes (Magma, Viridis, etc.), and font sizes
- **Export**: Download high-resolution plots

### 10. Gating Tab
- **Interactive Gating**: Select two features (Protein/Gene)
- **Lasso Select**: Draw around cells of interest to subset the Seurat object
- **Reset**: Revert to original dataset

### 11. VDJ Analysis Tab
- **Upload VDJ Data**: Load `filtered_contig_annotations.csv` (10x Genomics format)
- **Clonal Expansion**: Visualize clone size distribution across clusters
- **Diversity Analysis**: Compute and compare diversity indices (Shannon, Simpson)
- **Export**: Download plots and tables

### 12. Export Tab
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
  - Analysis: `UCell`, `DT`, `ggrepel`, `NMF`
  - VDJ Analysis: `scRepertoire`
  - Pathway Enrichment: `clusterProfiler`, `enrichplot`, `fgsea`, `msigdbr`, `DOSE`, `igraph`
  - Gene Annotation: `org.Hs.eg.db` (Human), `org.Mm.eg.db` (Mouse), `AnnotationDbi`
  - Utilities: `colourpicker`, `cowplot`, `plotly`

### For Docker
- Docker and Docker Compose

## Tips & Tricks

- **Large Datasets**: SCpubr automatically uses rasterization for >50,000 cells
- **Many Cell Types**: Manual color picker is scrollable for >20 levels
- **Feature Search**: Start typing in feature boxes for instant autocomplete
- **Split Plots**: Use `split.by` to facet plots by experimental conditions
- **Orientation**: Flip DotPlot/ViolinPlot to horizontal for better readability
- **NMF Analysis**: Start with k=5-10 factors; higher values find more granular programs but take longer
- **Pathway Enrichment**: Use ORA for discrete gene sets, GSEA for ranked gene lists
- **Gene Direction**: For ORA, filter by upregulated or downregulated genes to find direction-specific pathways
- **Heatmap Scaling**: Use row Z-score when comparing genes with different expression magnitudes
- **Ensembl IDs**: Convert Ensembl IDs to gene symbols before running analyses for better interpretability

## Troubleshooting

**App won't start**: Check that all dependencies are installed via `setup.R`

**Upload fails**: Ensure your file is a valid Seurat object and <10GB

**Plots not updating**: Click the plot to activate it and check its settings tab

**SCpubr errors**: Some metadata combinations may not be supported; try Standard style

## Project Structure

```
seurat_shiny_app/
├── app.R                    # Main application
├── enrichment_backend.R     # Pathway enrichment logic
├── setup.R                  # Dependency installer
├── Dockerfile               # Docker image definition
├── docker-compose.yml       # Docker orchestration
├── DOCKER_INSTRUCTIONS.md   # Docker deployment guide
└── README.md                # This file
```

## Citation

If you use this app in your research, please cite:
- **Seurat**: Hao et al., Cell (2021)
- **SCpubr**: Blanco-Carmona, E. (2022)
- **UCell**: Andreatta & Carmona, Nat Comput Sci (2021)
- **clusterProfiler**: Wu et al., The Innovation (2021); Yu et al., OMICS (2012)
- **NMF**: Gaujoux & Seoighe, BMC Bioinformatics (2010)
- **msigdbr**: Dolgalev, I. (2022)

## License

MIT License - see LICENSE file for details
