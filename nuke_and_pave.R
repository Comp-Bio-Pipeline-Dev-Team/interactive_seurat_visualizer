
# Get list of packages that are too new
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(version = "3.19", ask = FALSE, update = FALSE)
invalid <- BiocManager::valid(fix=FALSE)

if (!isTRUE(invalid) && !is.null(invalid$too_new)) {
  pkgs_to_remove <- rownames(invalid$too_new)
  print(paste("Removing", length(pkgs_to_remove), "conflicting packages..."))
  
  # Remove them
  remove.packages(pkgs_to_remove)
  print("Removal complete.")
}

print("Installing clusterProfiler (clean)...")
BiocManager::install("clusterProfiler", version="3.19", ask=FALSE, force=TRUE, update=TRUE)
