# Install CRAN packages
install.packages(c(
  "shiny", 
  "Seurat", 
  "ggplot2", 
  "patchwork", 
  "shinythemes", 
  "colourpicker", 
  "cowplot", 
  "plotly", 
  "shiny",
  "Seurat",
  "ggplot2",
  "patchwork",
  "shinythemes",
  "colourpicker",
  "cowplot",
  "plotly",
  "remotes",
  "BiocManager",
  "MetBrewer",
  "viridis",
  "RColorBrewer",
  "scales",
  "DT",
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
  "org.Mm.eg.db",     # Mouse annotations
  "DOSE"              # Disease ontology
# Install GitHub packages
remotes::install_github("enblacar/SCpubr")
