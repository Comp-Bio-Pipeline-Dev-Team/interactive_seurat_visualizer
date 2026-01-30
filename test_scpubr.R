library(Seurat)
library(ggplot2)
# Check if SCpubr involves
if (!requireNamespace("SCpubr", quietly = TRUE)) {
  cat("SCpubr not installed\n")
  quit(status=0)
}

obj <- readRDS("test_seurat.rds")
# Test standard ViolinPlot
tryCatch({
  p <- SCpubr::do_ViolinPlot(sample = obj, features = "IGLL5")
  cat("Basic SCpubr ViolinPlot: Success\n")
}, error = function(e) {
  cat("Basic SCpubr ViolinPlot: Failed - ", e$message, "\n")
})

# Test with split.plot = TRUE but no split.by
tryCatch({
  p <- SCpubr::do_ViolinPlot(sample = obj, features = "IGLL5", split.plot = TRUE)
  cat("Split=TRUE NoSplitBy: Success\n")
}, error = function(e) {
  cat("Split=TRUE NoSplitBy: Failed - ", e$message, "\n")
})
