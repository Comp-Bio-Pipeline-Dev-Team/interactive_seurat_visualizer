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
    p("Specify the data type for each metadata column. This affects how colors are applied in plots.", 
      style = "color: #7f8c8d; font-size: 14px;"),
    hr(),
    div(class = "metadata-table",
      p(strong("Note:"), "Metadata type editing will be implemented in the next iteration.", 
        style = "color: #95a5a6; font-style: italic;")
    )
  )
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
