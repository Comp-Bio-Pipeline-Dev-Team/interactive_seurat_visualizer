# Landing Page UI Module
# This module provides the upload interface shown before the main app

landing_page_ui <- function() {
  fluidPage(
    theme = shinytheme("flatly"),
    tags$head(tags$style(HTML("
      .landing-container {
        max-width: 900px;
        margin: 80px auto;
        padding: 40px;
        background: white;
        border-radius: 8px;
        box-shadow: 0 2px 15px rgba(0,0,0,0.1);
      }
      .upload-area {
        border: 3px dashed #3498db;
        border-radius: 8px;
        padding: 40px;
        text-align: center;
        background: #ecf0f1;
        margin-bottom: 30px;
      }
      .preview-box {
        background: #f8f9fa;
        border: 1px solid #dee2e6;
        border-radius: 5px;
        padding: 20px;
        margin: 20px 0;
      }
      .metadata-table {
        margin-top: 20px;
      }
      .proceed-btn {
        margin-top: 30px;
        text-align: center;
      }
    "))),
    
    div(class = "landing-container",
      h1("Seurat Interactive Visualizer", 
         style = "text-align: center; color: #2c3e50; margin-bottom: 10px;"),
      p("Upload your Seurat object to begin analysis", 
        style = "text-align: center; color: #7f8c8d; margin-bottom: 30px;"),
      hr(),
      
      # Upload Section
      div(class = "upload-area",
        icon("upload", style = "font-size: 48px; color: #3498db; margin-bottom: 15px;"),
        h3("Upload Your Seurat Object", style = "margin-bottom: 20px;"),
        fileInput("seurat_file", NULL, accept = ".rds", 
                  buttonLabel = "Browse Files", 
                  placeholder = "Select .rds file",
                  width = "100%"),
        p("Supported format: .rds (Seurat object)", 
          style = "color: #7f8c8d; font-size: 14px; margin-top: 10px;")
      ),
      
      # Object Preview (shown after upload)
      uiOutput("object_preview"),
      
      # Metadata Type Specification (shown after upload)
      uiOutput("metadata_type_ui"),
      
      # Ensembl Conversion (shown after upload)
      uiOutput("ensembl_conversion_ui"),
      
      # Proceed Button (shown after upload)
      uiOutput("proceed_button_ui")
    )
  )
}
