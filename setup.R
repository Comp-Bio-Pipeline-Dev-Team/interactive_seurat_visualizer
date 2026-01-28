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
  "remotes", 
  "BiocManager", 
  "MetBrewer", 
  "viridis",
  "RColorBrewer",
  "scales",
  "DT",
  "ggrepel"
))

# Install Bioconductor packages
BiocManager::install("UCell", update = FALSE)

# Install GitHub packages
remotes::install_github("enblacar/SCpubr")
