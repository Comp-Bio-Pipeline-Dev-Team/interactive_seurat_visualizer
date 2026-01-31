
# pseudobulk_backend.R
# Backend logic for Pseudobulk Heatmap

# Render Pseudobulk Heatmap
# Render Pseudobulk Heatmap
output$pseudobulk_heatmap <- renderPlot({
  req(data_object(), input$pb_group_by)
  
  obj <- data_object()
  group_by <- input$pb_group_by
  
  pb <- AggregateExpression(obj, group.by = group_by, return.seurat = FALSE)
  
  mat <- pb$RNA
  if (is.null(mat)) mat <- pb[[1]] 
  
  features <- VariableFeatures(obj)
  if (length(features) == 0) features <- rownames(mat)[1:50]
  features <- intersect(features, rownames(mat))
  
  mat_sub <- mat[features, , drop=FALSE]
  
  # Scaling
  if (input$pb_scaling == "row") {
    mat_scaled <- t(scale(t(mat_sub)))
    # Clip extreme values for better visualization
    mat_scaled[mat_scaled > 4] <- 4
    mat_scaled[mat_scaled < -4] <- -4
  } else {
    # For 'none', use log1p so it's viewable
    mat_scaled <- log1p(mat_sub)
  }
  
  # Palette Logic
  if (input$pb_palette == "Viridis") {
    color_fun <- viridis::viridis(100)
  } else if (input$pb_palette == "Magma") {
    color_fun <- viridis::magma(100)
  } else if (input$pb_palette == "RdBu") {
    color_fun <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))(100)
  } else {
    # Default RdYlBu
    color_fun <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100)
  }
  
  pheatmap::pheatmap(mat_scaled, 
                     main = paste("Pseudobulk Expression by", group_by),
                     show_rownames = length(features) < 50,
                     cluster_rows = input$pb_cluster_rows,
                     cluster_cols = input$pb_cluster_cols,
                     color = color_fun)
})

output$download_pb_heatmap <- downloadHandler(
  filename = function() { paste0("pseudobulk_heatmap_", Sys.Date(), ".pdf") },
  content = function(file) {
    req(data_object(), input$pb_group_by)
    # Regenerate Matrix (simplified for export)
    obj <- data_object()
    group_by <- input$pb_group_by
    pb <- AggregateExpression(obj, group.by = group_by, return.seurat = FALSE)
    mat <- pb$RNA
    if (is.null(mat)) mat <- pb[[1]]
    features <- VariableFeatures(obj)
    if (length(features) == 0) features <- rownames(mat)[1:50] 
    features <- intersect(features, rownames(mat))
    mat_sub <- mat[features, , drop=FALSE]
    if (input$pb_scaling == "row") mat_scaled <- t(scale(t(mat_sub))) else mat_scaled <- mat_sub
    
    color_fun <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
    if (input$pb_palette == "Viridis") color_fun <- viridis::viridis(100)
    if (input$pb_palette == "Magma") color_fun <- viridis::magma(100)
    
    pdf(file, width = input$pb_export_w, height = input$pb_export_h)
    if (requireNamespace("pheatmap", quietly=TRUE)) {
      pheatmap::pheatmap(mat_scaled, 
                         main = paste("Pseudobulk Expression by", group_by),
                         show_rownames = length(features) < 50,
                         cluster_rows = input$pb_cluster_rows,
                         cluster_cols = input$pb_cluster_cols,
                         color = color_fun)
    } else {
        heatmap(mat_scaled)
    }
    dev.off()
  }
)
