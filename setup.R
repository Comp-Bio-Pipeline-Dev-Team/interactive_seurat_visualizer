# Install CRAN packages
install.packages(c(
  "shiny", 
  "Seurat", 
  "ggplot2", 
  "patchwork", 
  "shinythemes", 
  "shinyjs",
  "colourpicker", 
  "cowplot", 
  "plotly", 
  "remotes", 
  "BiocManager", 
  "MetBrewer", 
  "pheatmap",
  "viridis",
  "RColorBrewer",
  "scales",
  "DT",
  "ggrepel",
  "msigdbr",  # MSigDB gene sets
  "igraph"    # For network plots
))

# Install Bioconductor packages
BiocManager::install(c(
  "UCell",
  "clusterProfiler",  # Pathway enrichment
  "enrichplot",       # Enrichment visualization
  "fgsea",            # Fast GSEA
  "org.Hs.eg.db",     # Human annotations
  "DOSE",             # Disease ontology
  "scRepertoire",     # VDJ Analysis
  "zellkonverter",    # H5AD Support
  "SingleCellExperiment",
  "NMF"               # Gene Programs
), update = FALSE)

# Install GitHub packages
remotes::install_github("enblacar/SCpubr")
