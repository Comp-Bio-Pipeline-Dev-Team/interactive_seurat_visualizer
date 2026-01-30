# Initialize renv for this project
# This script sets up renv and creates a lockfile for reproducible builds

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv (this will scan the project and create renv.lock)
renv::init()

# Snapshot the current state of packages
renv::snapshot()

cat("renv initialized successfully!\n")
cat("renv.lock file created with current package versions.\n")
