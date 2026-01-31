#!/bin/bash
echo "Checking for NMF package..."
R -e "if (!require('NMF', quietly=TRUE)) install.packages('NMF', repos='https://cloud.r-project.org')"

echo "Starting background environment restoration..."
nohup Rscript /app/restore_env_background.R > /app/restore.log 2>&1 &

echo "Starting Shiny App..."
exec R -e "shiny::runApp('/app', host='0.0.0.0', port=3838)"
