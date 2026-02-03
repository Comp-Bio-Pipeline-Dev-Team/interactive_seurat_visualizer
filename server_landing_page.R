# Landing Page Server Logic
# This file is sourced inside the server function in app.R
# All code here has access to input, output, session, and reactive values

# Object Preview UI
output$object_preview <- renderUI({
  req(seurat_obj())
  obj <- seurat_obj()
  
  div(class = "preview-box",
    h4("Object Preview", style = "color: #2c3e50;"),
    hr(),
    fluidRow(
      column(3,
        div(style = "text-align: center;",
          h2(ncol(obj), style = "color: #3498db; margin: 0;"),
          p("Cells", style = "color: #7f8c8d; margin: 0;")
        )
      ),
      column(3,
        div(style = "text-align: center;",
          h2(nrow(obj), style = "color: #e74c3c; margin: 0;"),
          p("Features", style = "color: #7f8c8d; margin: 0;")
        )
      ),
      column(3,
        div(style = "text-align: center;",
          h2(length(Assays(obj)), style = "color: #2ecc71; margin: 0;"),
          p("Assays", style = "color: #7f8c8d; margin: 0;")
        )
      ),
      column(3,
        div(style = "text-align: center;",
          h2(length(Reductions(obj)), style = "color: #f39c12; margin: 0;"),
          p("Reductions", style = "color: #7f8c8d; margin: 0;")
        )
      )
    ),
    hr(),
    p(strong("Assays:"), paste(Assays(obj), collapse = ", ")),
    p(strong("Reductions:"), paste(Reductions(obj), collapse = ", "))
  )
})

# Metadata Type Inference and Specification UI
output$metadata_type_ui <- renderUI({
  req(seurat_obj())
  
  div(class = "preview-box",
    h4("Metadata Column Types", style = "color: #2c3e50;"),
    p("Specify the data type for each metadata column. Categorical columns use discrete colors, numerical columns use continuous gradients.", 
      style = "color: #7f8c8d; font-size: 14px;"),
    hr(),
    div(class = "metadata-table",
      DT::dataTableOutput("metadata_type_table"),
      br(),
      actionButton("apply_metadata_types", "Apply Changes", 
                   icon = icon("check"), 
                   class = "btn-primary",
                   style = "margin-top: 10px;")
    )
  )
})

# Helper function to infer metadata types
infer_metadata_types <- function(obj) {
  meta <- obj@meta.data
  types <- sapply(colnames(meta), function(col) {
    values <- meta[[col]]
    
    # If already a factor, it's categorical
    if (is.factor(values)) return("Categorical")
    
    # If character, it's categorical
    if (is.character(values)) return("Categorical")
    
    # If numeric, check unique values
    if (is.numeric(values)) {
      n_unique <- length(unique(values))
      n_total <- length(values)
      
      # If few unique values relative to total, likely categorical
      if (n_unique < 20 && n_unique / n_total < 0.1) {
        return("Categorical")
      }
      return("Numerical")
    }
    
    # Default to categorical
    return("Categorical")
  })
  
  data.frame(
    Column = names(types),
    Type = as.character(types),
    stringsAsFactors = FALSE
  )
}

# Render metadata type table
output$metadata_type_table <- DT::renderDataTable({
  req(seurat_obj())
  
  type_df <- infer_metadata_types(seurat_obj())
  
  DT::datatable(
    type_df,
    selection = 'none',
    editable = list(target = 'cell', disable = list(columns = c(0))),
    options = list(
      pageLength = 10,
      dom = 't',
      ordering = FALSE
    ),
    rownames = FALSE
  )
})

# Apply metadata type changes
observeEvent(input$apply_metadata_types, {
  req(seurat_obj())
  
  # Get edited table data
  type_info <- input$metadata_type_table_cell_edit
  
  if (!is.null(type_info)) {
    obj <- seurat_obj()
    
    # Apply type conversions
    for (i in seq_len(nrow(type_info))) {
      row <- type_info$row[i]
      col_name <- infer_metadata_types(obj)$Column[row]
      new_type <- type_info$value[i]
      
      if (new_type == "Categorical") {
        # Convert to factor
        obj@meta.data[[col_name]] <- as.factor(obj@meta.data[[col_name]])
      } else if (new_type == "Numerical") {
        # Convert to numeric (if possible)
        obj@meta.data[[col_name]] <- as.numeric(as.character(obj@meta.data[[col_name]]))
      }
    }
    
    seurat_obj(obj)
    showNotification("Metadata types updated successfully!", type = "message", duration = 3)
  } else {
    showNotification("No changes detected", type = "warning", duration = 2)
  }
})


# Ensembl Conversion UI
output$ensembl_conversion_ui <- renderUI({
  req(seurat_obj())
  
  div(class = "preview-box",
    h4("Gene ID Conversion (Optional)", style = "color: #2c3e50;"),
    p("Convert Ensembl IDs to gene symbols for easier interpretation.", 
      style = "color: #7f8c8d; font-size: 14px;"),
    hr(),
    radioButtons("ensembl_species", "Species:", 
                 choices = c("Human", "Mouse"), 
                 inline = TRUE,
                 selected = "Human"),
    actionButton("convert_ensembl_landing", "Convert Ensembl IDs to Symbols", 
                 icon = icon("dna"), 
                 class = "btn-info")
  )
})

# Handle Ensembl ID conversion
observeEvent(input$convert_ensembl_landing, {
  req(seurat_obj())
  
  species <- input$ensembl_species
  db_pkg <- if(species == "Human") "org.Hs.eg.db" else "org.Mm.eg.db"
  
  # Check if required packages are installed
  if (!requireNamespace(db_pkg, quietly = TRUE)) {
    showNotification(paste(db_pkg, "is not installed. Please install it to use this feature."), 
                     type="error", duration = 5)
    return()
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    showNotification("AnnotationDbi is not installed. Please install it to use this feature.", 
                     type="error", duration = 5)
    return()
  }
  
  withProgress(message = "Converting Ensembl IDs...", value = 0, {
    tryCatch({
      # Load the package dynamically
      library(db_pkg, character.only = TRUE)
      org_db <- get(db_pkg)
      
      obj <- seurat_obj()
      current_ids <- rownames(obj)
      
      # Detect if IDs are Ensembl format (start with ENS)
      ensembl_pattern <- grepl("^ENS", current_ids)
      if (sum(ensembl_pattern) == 0) {
        showNotification("No Ensembl IDs detected in feature names.", 
                        type="warning", duration = 4)
        return()
      }
      
      # Map IDs
      incProgress(0.3, detail = paste("Mapping", species, "Ensembl to Symbols"))
      mapped_symbols <- AnnotationDbi::mapIds(
        org_db,
        keys = current_ids,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
      
      # Identify features that successfully mapped
      mapped_indices <- !is.na(mapped_symbols)
      num_mapped <- sum(mapped_indices)
      
      if (num_mapped == 0) {
        showNotification("No Ensembl IDs could be mapped to gene symbols.", 
                        type="warning", duration = 4)
        return()
      }
      
      # Create new names: Use symbol if mapped, else keep original ID
      incProgress(0.6, detail = "Creating new feature names")
      new_names <- current_ids
      new_names[mapped_indices] <- mapped_symbols[mapped_indices]
      
      # Handle duplicates by making names unique
      new_names <- make.unique(new_names)
      
      # Update feature names in all assays
      incProgress(0.8, detail = "Updating Seurat object")
      for (assay_name in Assays(obj)) {
        assay_obj <- obj[[assay_name]]
        rownames(assay_obj) <- new_names
        obj[[assay_name]] <- assay_obj
      }
      
      # Update seurat_obj
      seurat_obj(obj)
      
      showNotification(
        paste("Successfully converted", num_mapped, "Ensembl IDs to gene symbols!"),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      showNotification(paste("Error during conversion:", e$message), 
                      type="error", duration = 5)
    })
  })
})


# Proceed Button UI
output$proceed_button_ui <- renderUI({
  req(seurat_obj())
  
  div(class = "proceed-btn",
    actionButton("proceed_to_app", "Proceed to App", 
                 icon = icon("arrow-right"),
                 class = "btn-success btn-lg",
                 style = "padding: 15px 40px; font-size: 18px;")
  )
})

# Handle Proceed to App button
observeEvent(input$proceed_to_app, {
  req(seurat_obj())
  
  # Set app_ready to TRUE to show main app
  app_ready(TRUE)
})
