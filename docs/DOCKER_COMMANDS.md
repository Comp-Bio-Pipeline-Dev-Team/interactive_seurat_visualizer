# Docker Commands (without docker-compose)

## Building the Image

```bash
cd /home/tonya/Downloads/github/interactive_seurat_visualizer
sudo docker build -t seurat-shiny-app .
```

## Running the Container

### Basic run:
```bash
sudo docker run -p 3838:3838 seurat-shiny-app
```

### Run in background (detached):
```bash
sudo docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app
```

### Run with data volume mounted:
```bash
sudo docker run -d -p 3838:3838 \
  -v $(pwd)/data:/data \
  --name seurat-app \
  seurat-shiny-app
```

## Accessing the App

Open browser to: **http://localhost:3838**

## Managing the Container

### View logs:
```bash
sudo docker logs seurat-app
```

### Follow logs in real-time:
```bash
sudo docker logs -f seurat-app
```

### Stop the container:
```bash
sudo docker stop seurat-app
```

### Start stopped container:
```bash
sudo docker start seurat-app
```

### Remove container:
```bash
sudo docker rm seurat-app
```

### Remove container and image:
```bash
sudo docker stop seurat-app
sudo docker rm seurat-app
sudo docker rmi seurat-shiny-app
```

## Rebuilding After Changes

```bash
# Stop and remove old container
sudo docker stop seurat-app
sudo docker rm seurat-app

# Rebuild image
sudo docker build -t seurat-shiny-app .

# Run new container
sudo docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app
```

## Quick Rebuild Script

Save this as `rebuild.sh`:
```bash
#!/bin/bash
sudo docker stop seurat-app 2>/dev/null
sudo docker rm seurat-app 2>/dev/null
sudo docker build -t seurat-shiny-app .
sudo docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app
echo "App running at http://localhost:3838"
```

Make executable and run:
```bash
chmod +x rebuild.sh
./rebuild.sh
```
