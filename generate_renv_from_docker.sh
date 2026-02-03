#!/bin/bash
# Script to generate renv.lock from running Docker container

set -e

CONTAINER_NAME="seurat-app"
OUTPUT_FILE="renv.lock"

echo "Generating renv.lock from Docker container: $CONTAINER_NAME"

# Step 1: Execute renv::snapshot() inside the container
echo "Step 1: Running renv::snapshot() in container..."
docker exec $CONTAINER_NAME R -e "
  # Initialize renv if needed
  if (!file.exists('renv.lock')) {
    renv::init(bare = TRUE)
  }
  
  # Snapshot all installed packages
  renv::snapshot(prompt = FALSE)
  
  cat('renv.lock generated successfully\n')
"

# Step 2: Copy the generated renv.lock from container to host
echo "Step 2: Copying renv.lock from container to host..."
docker cp $CONTAINER_NAME:/app/renv.lock ./$OUTPUT_FILE

# Step 3: Fix ownership
echo "Step 3: Fixing file ownership..."
sudo chown $USER:$USER ./$OUTPUT_FILE

echo "✅ Done! renv.lock has been generated and saved to: $OUTPUT_FILE"
echo ""
echo "File stats:"
wc -l $OUTPUT_FILE
ls -lh $OUTPUT_FILE
