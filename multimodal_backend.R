
# multimodal_backend.R
# Backend logic for Multimodal (ADT/CITE-seq) Analysis

# Reactive to check for ADT assay
has_adt <- reactive({
  req(data_object())
  "ADT" %in% Assays(data_object()) || "Protein" %in% Assays(data_object())
})

adt_assay_name <- reactive({
  req(data_object())
  if ("ADT" %in% Assays(data_object())) return("ADT")
  if ("Protein" %in% Assays(data_object())) return("Protein")
  return(NULL)
})

# Update UI based on ADT presence
observe({
  req(has_adt())
  assay <- adt_assay_name()
  
  # Update feature choices for ADT
  features <- rownames(data_object()[[assay]])
  updateSelectizeInput(session, "adt_feature", choices = features, server = TRUE)
  updateSelectizeInput(session, "adt_feature_x", choices = features, server = TRUE)
  
  # Update RNA features for co-expression (using variable features to limit list size)
  rna_features <- VariableFeatures(data_object())
  if (length(rna_features) == 0) rna_features <- rownames(data_object())[1:100] # Fallback
  updateSelectizeInput(session, "rna_feature_y", choices = rna_features, server = TRUE)
})

# Render ADT FeaturePlot
output$adt_feature_plot <- renderPlot({
  req(data_object(), input$adt_feature)
  
  obj <- data_object()
  feat <- input$adt_feature
  assay <- adt_assay_name()
  
  # Ensure we are pulling from the correct slot
  DefaultAssay(obj) <- assay
  
  p <- FeaturePlot(obj, features = feat, pt.size = input$pt_size_adt, order = TRUE) +
       scale_color_viridis_c(option = input$adt_palette) +
       theme_cowplot() +
       ggtitle(paste("Protein Expression:", feat))
  
  return(p)
})


# Render Co-expression Plot (RNA vs ADT)
output$adt_rna_coexpression <- renderPlot({
  req(data_object(), input$adt_feature_x, input$rna_feature_y)
  
  obj <- data_object()
  prot <- input$adt_feature_x
  gene <- input$rna_feature_y
  assay <- adt_assay_name()
  
  # Fetch data
  # Check if protein and gene are same to avoid error? Usually they differ.
  prot_vals <- FetchData(obj, vars = paste0(assay, "_", prot))
  gene_vals <- FetchData(obj, vars = paste0("rna_", gene))
  
  df <- data.frame(Protein = prot_vals[,1], Gene = gene_vals[,1])
  
  p <- ggplot(df, aes(x = Protein, y = Gene)) +
       geom_point(alpha = 0.6, color = "#2c3e50", size = 0.5) +
       geom_smooth(method = "lm", color = "#e74c3c", se = FALSE) +
       theme_bw() +
       labs(x = paste("Protein:", prot), y = paste("Gene:", gene)) +
       ggtitle(paste("Co-expression:", prot, "vs", gene))
       
  return(p)
})

# Download Handlers
output$download_adt_feature <- downloadHandler(
  filename = function() { paste0("adt_feature_", input$adt_feature, ".", input$adt_export_format) },
  content = function(file) {
    req(data_object(), input$adt_feature)
    obj <- data_object()
    DefaultAssay(obj) <- adt_assay_name()
    p <- FeaturePlot(obj, features = input$adt_feature, pt.size = input$pt_size_adt, order=TRUE) + 
         scale_color_viridis_c(option = input$adt_palette)
    ggsave(file, plot = p, device = input$adt_export_format, width = input$adt_export_width, height = input$adt_export_height)
  }
)

output$download_adt_coexp <- downloadHandler(
  filename = function() { paste0("adt_coexp_", input$adt_feature_x, "_", input$rna_feature_y, ".", input$adt_export_format) },
  content = function(file) {
    req(data_object(), input$adt_feature_x, input$rna_feature_y)
    # Re-generate plot (ideal to refactor into reactive, but copy-paste safe for now)
    obj <- data_object()
    prot <- input$adt_feature_x
    gene <- input$rna_feature_y
    assay <- adt_assay_name()
    prot_vals <- FetchData(obj, vars = paste0(assay, "_", prot))
    gene_vals <- FetchData(obj, vars = paste0("rna_", gene))
    df <- data.frame(Protein = prot_vals[,1], Gene = gene_vals[,1])
    p <- ggplot(df, aes(x = Protein, y = Gene)) +
         geom_point(alpha = 0.6, color = "#2c3e50", size = 0.5) +
         geom_smooth(method = "lm", color = "#e74c3c", se = FALSE) +
         theme_bw() +
         labs(x = paste("Protein:", prot), y = paste("Gene:", gene)) +
         ggtitle(paste("Co-expression:", prot, "vs", gene))
    
    ggsave(file, plot = p, device = input$adt_export_format, width = input$adt_export_width, height = input$adt_export_height)
  }
)
