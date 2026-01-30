FROM rocker/r-ver:4.4.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    libgit2-dev \
    libhdf5-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    pandoc \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Install R packages - Core & CRAN
# Split installation to leverage caching layers
RUN R -e "install.packages(c('shiny', 'shinythemes', 'ggplot2', 'patchwork', 'cowplot', 'plotly', 'colourpicker', 'MetBrewer', 'viridis', 'RColorBrewer', 'remotes', 'BiocManager', 'scales', 'devtools', 'DT', 'ggrepel'), repos='https://cloud.r-project.org')"

# Install Seurat (separate step closely following core deps)
RUN R -e "install.packages('Seurat', repos='https://cloud.r-project.org')"

# Install UCell from Bioconductor
RUN R -e "BiocManager::install('UCell', update=FALSE)"

# Install SCpubr from GitHub
RUN R -e "remotes::install_github('enblacar/SCpubr')"

# Copy application files
COPY app.R /app/
COPY color_utils.R /app/
COPY plot_cluster_distribution.R /app/
COPY plot_dimension_reduction.R /app/
COPY plot_feature.R /app/
COPY plot_violin.R /app/
COPY plot_dot.R /app/
COPY setup.R /app/
COPY generate_dummy_seurat.R /app/
COPY README.md /app/

# Expose port for Shiny
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
