
# vdj_backend.R
# Backend logic for VDJ/Repertoire Analysis using scRepertoire

# Reactive value for VDJ data
vdj_data <- reactiveVal(NULL)

# Helper: Check if scRepertoire is installed
has_screp <- requireNamespace("scRepertoire", quietly = TRUE)

# Handle file upload
observeEvent(input$vdj_upload, {
  req(input$vdj_upload)
  
  if (!has_screp) {
    showNotification("scRepertoire package is missing. Please install it.", type = "error")
    return()
  }
  
  tryCatch({
    # Read the text file
    # SCRepertoire expects a list of contig dataframes for combineTCR/BCR
    # Here we assume a single sample upload for simplicity, or a merged file
    df <- read.csv(input$vdj_upload$datapath)
    
    # Check minimal columns (barcode, chain, v_gene, j_gene, cdr3)
    if (!all(c("barcode", "chain", "cdr3") %in% colnames(df))) {
      showNotification("Invalid contig file. Must have barcode, chain, cdr3 columns.", type = "error")
      return()
    }
    
    # Process using scRepertoire
    # We guess ID based on Seurat object if possible
    combined <- scRepertoire::combineTCR(df, samples = "Sample1", ID = "ID1")
    
    vdj_data(combined)
    showNotification("VDJ data loaded successfully.", type = "message")
    
    # Attempt integration with Seurat object
    if (!is.null(data_object())) {
      obj <- data_object()
      # scRepertoire::combineExpression adds metadata to the object
      # Note: match by barcodes is critical.
      obj <- scRepertoire::combineExpression(combined, obj, cloneCall="gene+nt")
      data_object(obj) # Update the main object
      showNotification("VDJ data integrated with Seurat object.", type = "message")
    }
    
  }, error = function(e) {
    showNotification(paste("Error loading VDJ:", e$message), type = "error")
  })
})

# Render Clonal Overlay
output$vdj_overlay <- renderPlot({
  req(data_object())
  if (!has_screp) return(NULL)
  
  # Check if clone info is in meta.data
  if (!"cloneType" %in% colnames(data_object()@meta.data)) {
    return(NULL)
  }
  
  scRepertoire::clonalOverlay(data_object(), reduction = "umap", 
                              freq.cutpoints = 1, bins = 10, facet = "chain") +
    theme_cowplot()
})


# Render V Gene Usage
output$vdj_vgene <- renderPlot({
  req(vdj_data())
  if (!has_screp) return(NULL)
  
  scRepertoire::vizGenes(vdj_data(), gene = "V", chain = "TRB", plot = "heatmap", scale = TRUE)
})

# Update Group By choices for diversity
observe({
    req(data_object())
    # Update group by choices for diversity
    numeric_cols <- sapply(data_object()@meta.data, is.numeric)
    choices <- colnames(data_object()@meta.data)
    updateSelectInput(session, "vdj_group_by", choices = choices, selected = "orig.ident")
})

# Render Diversity Plot
output$vdj_diversity_plot <- renderPlot({
    req(data_object(), input$vdj_group_by, input$vdj_diversity_metric)
    if (!has_screp) return(NULL)
    
    # Check if necessary columns exist (done by scRepertoire, but good to ensure object is set)
    obj <- data_object()
    
    # scRepertoire::clonalDiversity
    # Note: Requires the index to be valid
    scRepertoire::clonalDiversity(obj, 
                                  cloneCall = "gene+nt", 
                                  group.by = input$vdj_group_by, 
                                  n.boots = input$vdj_boot_n, 
                                  metrics = input$vdj_diversity_metric) +
        theme_cowplot() +
        ggtitle(paste("Clonal Diversity:", input$vdj_diversity_metric))
})

# Download Handler for VDJ
output$download_vdj_plot <- downloadHandler(
  filename = function() { paste0("vdj_plot_", Sys.Date(), ".png") },
  content = function(file) {
    req(data_object())
    if (!has_screp) return()
    
    p <- NULL
    if (input$vdj_viz_type == "Clonal Overlay") {
       if ("cloneType" %in% colnames(data_object()@meta.data)) {
         p <- scRepertoire::clonalOverlay(data_object(), reduction = "umap", 
                                  freq.cutpoints = 1, bins = 10, facet = "chain") +
              theme_cowplot()
       }
    } else if (input$vdj_viz_type == "V Gene Usage") {
       p <- scRepertoire::vizGenes(vdj_data(), gene = "V", chain = "TRB", plot = "heatmap", scale = TRUE)
    } else if (input$vdj_viz_type == "Clonal Diversity") {
       p <- scRepertoire::clonalDiversity(data_object(), 
                                  cloneCall = "gene+nt", 
                                  group.by = input$vdj_group_by, 
                                  n.boots = input$vdj_boot_n, 
                                  metrics = input$vdj_diversity_metric) +
        theme_cowplot() +
        ggtitle(paste("Clonal Diversity:", input$vdj_diversity_metric))
    }
    
    if (!is.null(p)) {
      ggsave(file, plot = p, width = 12, height = 8)
    }
  }
)
