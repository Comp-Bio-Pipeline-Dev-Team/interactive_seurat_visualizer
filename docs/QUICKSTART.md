# Quick Start Guide

## ✅ Setup Complete!

Your renv lockfile has been created successfully. The app is now ready to build with reproducible package versions.

## Build and Run (Using Docker Only)

### 1. Build the Docker image:
```bash
cd /home/tonya/Downloads/github/interactive_seurat_visualizer
sudo docker build -t seurat-shiny-app .
```

### 2. Run the container:
```bash
sudo docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app
```

### 3. Access the app:
Open browser to: **http://localhost:3838**

### 4. View logs:
```bash
sudo docker logs -f seurat-app
```

### 5. Stop the app:
```bash
sudo docker stop seurat-app
```

## Testing the ClusterDistrBar Fix

1. Upload a Seurat object (or use `test_seurat.rds`)
2. Select **Plot 1** → **Plot Type: ClusterDistrBar**
3. Choose X Axis and Fill variables
4. Set **Color Source: Manual**
5. Change colors using the pickers
6. ✅ **Verify:** Colors should apply immediately!

## Rebuilding After Code Changes

```bash
sudo docker stop seurat-app
sudo docker rm seurat-app
sudo docker build -t seurat-shiny-app .
sudo docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app
```

## Documentation

- **Docker Commands:** [DOCKER_COMMANDS.md](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/DOCKER_COMMANDS.md)
- **Renv Setup:** [RENV_SETUP.md](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/RENV_SETUP.md)
- **Refactoring Details:** [docs/](file:///home/tonya/Downloads/github/interactive_seurat_visualizer/docs/)

## Git Status

All changes committed to branch `dev_tb`:
- `75b79f4` - feat: Add unified color utilities module
- `da5801a` - docs: Add refactoring documentation to project
- Latest - build: Initialize renv with lockfile

Ready to push to remote!
