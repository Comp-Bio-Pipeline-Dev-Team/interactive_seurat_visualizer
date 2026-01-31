
# smart_nuke.R
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")

print("Checking for too-new packages WITHOUT triggering install error...")

# We need to manually remove packages that are "too new" for 3.19
# BiocManager::valid() checks against the *current* version. 
# We want to check against 3.19.
# But we can't switch to 3.19 without erroring.

# Strategy: 
# 1. Force remove known problematic packages (Bioconductor core)
# 2. Force install BiocVersion 3.19 (this defines the version)
# 3. Then run check and clean

core_pkgs <- c("BiocVersion", "BiocManager", "clusterProfiler", "enrichplot", "ggtree", "treeio", "GOSemSim", "DOSE", "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi", "Biobase", "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges")

print("Removing core Bioc packages...")
invisible(sapply(core_pkgs, function(x) try(remove.packages(x), silent=TRUE)))

print("Manually installing BiocVersion 3.19...")
install.packages("BiocVersion", repos="https://bioconductor.org/packages/3.19/bioc")

print("Now using BiocManager to heal...")
# This should now see we are on 3.19 (via BiocVersion) and downgrade others 
BiocManager::install(version = "3.19", ask = FALSE, checkBuilt = TRUE, update = TRUE, force = TRUE)

print("Installing clusterProfiler...")
BiocManager::install("clusterProfiler", version="3.19", ask=FALSE, force=TRUE)
