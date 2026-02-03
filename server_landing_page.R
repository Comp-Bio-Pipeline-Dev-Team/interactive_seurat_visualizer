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
                 class = "btn-info"),
    p(strong("Note:"), "Ensembl conversion will be implemented in the next iteration.", 
      style = "color: #95a5a6; font-style: italic; margin-top: 10px;")
  )
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
