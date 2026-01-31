
# qc_backend.R
# Backend logic for Exploratory QC

# Render QC Violin Plots
output$qc_violin <- renderPlot({
  req(data_object())
  
  # Check for percent.mt
  obj <- data_object()
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    # Attempt to calculate if missing
    if (any(grepl("^MT-", rownames(obj)))) {
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
        data_object(obj) # Update reactive object
    }
  }
  
  features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  features <- intersect(features, colnames(obj@meta.data))
  
  if (length(features) == 0) return(NULL)
  
  VlnPlot(obj, features = features, group.by = "orig.ident", pt.size = 0.1, ncol = length(features)) +
    theme_cowplot()
})

# Render QC Scatter Plot
output$qc_scatter <- renderPlot({
  req(data_object())
  obj <- data_object()
  
  feat1 <- "nCount_RNA"
  feat2 <- "percent.mt"
  
  if (all(c(feat1, feat2) %in% colnames(obj@meta.data))) {
      FeatureScatter(obj, feature1 = feat1, feature2 = feat2, group.by = "orig.ident") +
      theme_cowplot()
  } else {
      return(NULL)
  }
})
