# Renv Setup Instructions

## What is renv?

`renv` is an R package that creates reproducible environments by tracking exact package versions. This ensures that the Docker container builds with the same package  every time.

## Initial Setup (One-time)

Run this **once** to create the `renv.lock` file:

```bash
cd /path/to/interactive_seurat_visualizer

# Option 1: Using Docker (recommended)
docker run --rm -v $(pwd):/app -w /app rocker/r-ver:4.4.0 R -e "install.packages('renv'); renv::init(bare = TRUE); renv::snapshot()"

# Option 2: Using local R (if you have R installed)
Rscript init_renv.R
```

This will create:
- `renv.lock` - Lockfile with exact package versions
- `.Rprofile` - Activates renv when R starts
- `renv/` directory - renv infrastructure

## Using the renv-enabled Dockerfile

Once `renv.lock` exists, use the renv-enabled Dockerfile:

```bash
# Rename the current Dockerfile
mv Dockerfile Dockerfile.old

# Use the renv version
mv Dockerfile.renv Dockerfile

# Build with renv
docker-compose build
```

## Benefits

✅ **Reproducible builds** - Same package versions every time
✅ **Faster rebuilds** - Docker caches the package layer
✅ **Version control** - Track package changes in git
✅ **Conflict resolution** - Avoid "works on my machine" issues

## Updating Packages

When you add new packages to the app:

```r
# In R session
renv::snapshot()  # Update renv.lock with new packages
```

Then rebuild the Docker image to use the new versions.

## Troubleshooting

If renv restore fails in Docker:
```bash
# Clear Docker cache and rebuild
docker-compose build --no-cache
```
