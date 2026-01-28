# Script to generate a small dummy Seurat object
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat is not installed. Please install it first.")
}
library(Seurat)

# Use pbmc_small dataset which comes with Seurat
data("pbmc_small")

# Save it as RDS
saveRDS(pbmc_small, "test_seurat.rds")
print("Verified: test_seurat.rds created successfully.")
