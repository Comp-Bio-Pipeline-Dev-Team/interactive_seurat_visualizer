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
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Install renv
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"

# Copy renv.lock file
COPY renv.lock /app/

# Override Posit Package Manager (p3m.dev) with standard CRAN mirror.
# The renv.lock was created with p3m.dev URLs which produce malformed paths
# in some renv versions during Docker builds. cloud.r-project.org is stable.
ENV RENV_CONFIG_REPOS_OVERRIDE="https://cloud.r-project.org"

# Restore packages from lockfile using the CRAN mirror override
RUN R -e "renv::restore(prompt = FALSE)"

# Copy application files
COPY app.R /app/
COPY color_utils.R /app/
COPY plot_cluster_distribution.R /app/
COPY plot_dimension_reduction.R /app/
COPY plot_feature.R /app/
COPY plot_violin.R /app/
COPY plot_dot.R /app/
COPY ui_landing_page.R /app/
COPY server_landing_page.R /app/
COPY setup.R /app/
COPY generate_dummy_seurat.R /app/
COPY README.md /app/

# Expose port for Shiny
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
