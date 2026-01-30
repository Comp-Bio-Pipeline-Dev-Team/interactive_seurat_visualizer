# Docker Deployment Instructions

## Building the Docker Image

```bash
cd /Users/laurentgapin/.gemini/antigravity/scratch/seurat_shiny_app
docker build -t seurat-shiny-app .
```

## Running with Docker

### Option 1: Using docker run
```bash
docker run -p 3838:3838 seurat-shiny-app
```

Then open your browser to: `http://localhost:3838`

### Option 2: Using docker-compose (Recommended)
```bash
docker-compose up -d
```

To stop:
```bash
docker-compose down
```

## Accessing the App

Once running, open your browser to:
- **Local**: http://localhost:3838
- **Network**: http://YOUR_IP_ADDRESS:3838

## Volume Mounting for Data

To mount a directory with your Seurat objects:

```bash
docker run -p 3838:3838 -v /path/to/your/data:/data seurat-shiny-app
```

## Version Control

This Docker image includes:
- R 4.4.0
- Shiny and all visualization packages
- Seurat (latest from CRAN)
- UCell (from Bioconductor)
- All dependencies pinned for reproducibility

## Troubleshooting

### Check logs
```bash
docker logs <container-id>
```

### Interactive shell
```bash
docker run -it seurat-shiny-app /bin/bash
```

## Updating the App

1. Make changes to `app.R`
2. Rebuild the image:
   ```bash
   docker build -t seurat-shiny-app .
   ```
3. Restart the container:
   ```bash
   docker-compose restart
   ```
