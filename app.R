library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)
library(shinythemes)
library(colourpicker)
library(cowplot)
library(plotly)
library(MetBrewer)
library(viridis)
library(RColorBrewer)
library(DT)
library(ggrepel)

# Optional Libraries
has_ucell <- requireNamespace("UCell", quietly = TRUE)
has_scpubr <- requireNamespace("SCpubr", quietly = TRUE)
has_enrichment <- requireNamespace("clusterProfiler", quietly = TRUE) && 
                  requireNamespace("enrichplot", quietly = TRUE)
if (has_ucell) library(UCell)
if (has_scpubr) library(SCpubr)
if (has_enrichment) {
  library(clusterProfiler)
  library(enrichplot)
  library(fgsea)
  library(msigdbr)
  library(DOSE)
  library(igraph)
}

options(shiny.maxRequestSize = 10000 * 1024^2)

# --- Load Utility Modules ---
source("color_utils.R")
source("plot_cluster_distribution.R")

# --- UI Helper ---
plotControlUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4(paste("Settings for", id)), 
    selectInput(ns("plot_type"), "Plot Type", 
                choices = c("DimPlot", "FeaturePlot", "ViolinPlot", "DotPlot", "ClusterDistrBar")),
    
    uiOutput(ns("dynamic_ui")), 
    
    hr(),
    h5("Appearance & Colors"),
    fluidRow(
      column(6, selectInput(ns("title_align"), "Title Align", choices = c("Left"=0, "Center"=0.5, "Right"=1), selected=0.5)),
      column(6, uiOutput(ns("pt_size_ui")))
    ),
    
    # Orientation toggle for specific plot types
    uiOutput(ns("orientation_ui")),
    
    selectInput(ns("color_source"), "Color Source", choices = c("Default", "Palette", "Manual")),
    conditionalPanel(
      condition = "input.color_source == 'Palette'", ns = ns,
      selectInput(ns("palette_name"), "Palette", 
                  choices = list(
                    "Viridis" = c("viridis", "magma", "plasma", "inferno", "cividis"),
                    "RColorBrewer" = c("Set1", "Set2", "Set3", "Dark2", "Paired", "Pastel1", "Accent"),
                    "MetBrewer" = names(MetBrewer::MetPalettes)
                  ))
    ),
    conditionalPanel(
      condition = "input.color_source == 'Manual'", ns = ns,
      uiOutput(ns("manual_color_ui"))
    ),
    
    hr(),
    h5("Style & Advanced"),
    uiOutput(ns("plot_style_ui")),  # Dynamic based on plot type
    
    # SCpubr-specific controls
    conditionalPanel(
      condition = "input.plot_style == 'SCpubr'", ns = ns,
      numericInput(ns("scpubr_font_size"), "Base Font Size", value = 14, min = 8, max = 24),
      checkboxInput(ns("scpubr_border"), "Add Cell Borders", value = TRUE),
      selectInput(ns("scpubr_legend_pos"), "Legend Position", 
                  choices = c("bottom", "top", "left", "right", "none"), selected = "bottom"),
      numericInput(ns("scpubr_raster_dpi"), "Raster DPI (for large datasets)", value = 2048, min = 72, max = 4096),
      selectInput(ns("scpubr_viridis_dir"), "Viridis Direction", 
                  choices = c("Dark = Low (1)" = 1, "Dark = High (-1)" = -1), selected = 1)
    ),
    
    textInput(ns("custom_title"), "Custom Title", ""),
    numericInput(ns("title_size"), "Title Size", value = 16, min=8, max=30),
    checkboxInput(ns("show_legend"), "Show Legend", value = TRUE)
  )
}

ui <- navbarPage(
  title = "Seurat Interactive Visualizer",
  theme = shinytheme("flatly"),
  id = "main_nav",
  
  header = tags$head(tags$style(HTML("
    .plot-container { border: 1px solid #ddd; padding: 10px; margin: 10px; background: white; border-radius: 4px; }
    .active-plot { border: 3px solid #2c3e50 !important; box-shadow: 0 0 10px rgba(44,62,80,0.3); }
    .nav-tabs { font-weight: bold; }
  "))),
  
  tabPanel("Visualization",
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("seurat_file", "Upload Seurat (.rds or .h5ad)", accept = c(".rds", ".h5ad")),
      radioButtons("ensembl_species", "Species", choices = c("Human", "Mouse"), inline = TRUE),
      actionButton("convert_ensembl", "Convert Ensembl IDs to Symbols", icon = icon("dna"), class = "btn-info btn-block"),
      br(),
      
      # JavaScript to highlight active plot and switch tabs
      tags$script(HTML("
        $(document).on('shiny:inputchanged', function(event) {
          if (event.name === 'settings_tabs') {
            // Remove active class from all
            $('.plot-container').removeClass('active-plot');
            // Add to corresponding plot
            if (event.value === 'Plot 1') $('#plot1-container').addClass('active-plot');
            if (event.value === 'Plot 2') $('#plot2-container').addClass('active-plot');
            if (event.value === 'Plot 3') $('#plot3-container').addClass('active-plot');
            if (event.value === 'Plot 4') $('#plot4-container').addClass('active-plot');
          }
        });
        
        // Click on plot to activate and switch tab
        $(document).on('click', '.plot-container', function() {
          var plotId = $(this).attr('id');
          var tabName = '';
          if (plotId === 'plot1-container') tabName = 'Plot 1';
          if (plotId === 'plot2-container') tabName = 'Plot 2';
          if (plotId === 'plot3-container') tabName = 'Plot 3';
          if (plotId === 'plot4-container') tabName = 'Plot 4';
          
          if (tabName) {
            // Switch to the corresponding tab
            $('a[data-value=\"' + tabName + '\"]').tab('show');
            // Highlight the plot
            $('.plot-container').removeClass('active-plot');
            $(this).addClass('active-plot');
          }
        });
      ")),
      
      tabsetPanel(id = "settings_tabs",
        
        tabPanel("Subsetting",
          br(),
          h5("Preview (All Cells)"),
          plotOutput("preview_plot", height="250px"),
          hr(),
          selectInput("subset_col", "Meta Col", choices = NULL),
          selectizeInput("subset_levels", "Keep Levels", choices = NULL, multiple=TRUE),
          actionButton("apply_subset_meta", "Filter", class="btn-warning btn-sm"),
          hr(),
          actionButton("reset_data", "Reset", class="btn-danger")
        ),
        
        tabPanel("Signatures",
           br(),
           textAreaInput("sig_genes", "Genes (one per line)", rows=5, placeholder="CD3D\nCD3E\nCD8A"),
           textInput("sig_name", "Signature Name", "MySig"),
           actionButton("calc_ucell", "Calculate UCell", class="btn-success")
        ),
        
        tabPanel("Plot 1", plotControlUI("p1")),
        tabPanel("Plot 2", plotControlUI("p2")),
        tabPanel("Plot 3", plotControlUI("p3")),
        tabPanel("Plot 4", plotControlUI("p4")),
        
        tabPanel("Export",
           br(),
           selectInput("export_target", "Target", choices = c("All Plots", "Plot 1", "Plot 2", "Plot 3", "Plot 4")),
           selectInput("export_format", "Format", choices = c("png", "pdf", "jpg")),
           numericInput("export_w", "Width", 12), 
           numericInput("export_h", "Height", 10),
           downloadButton("download_plot", "Download")
        )
      )
    ),
    
    mainPanel(
      width = 9,
      textOutput("obj_info"),
      fluidRow(
        column(6, 
          div(style="text-align: center; font-weight: bold; color: #2c3e50; margin-bottom: 5px;", "Plot 1"),
          div(id="plot1-container", class="plot-container", plotOutput("plot1", height="400px"))
        ),
        column(6, 
          div(style="text-align: center; font-weight: bold; color: #2c3e50; margin-bottom: 5px;", "Plot 2"),
          div(id="plot2-container", class="plot-container", plotOutput("plot2", height="400px"))
        )
      ),
      fluidRow(
        column(6, 
          div(style="text-align: center; font-weight: bold; color: #2c3e50; margin-bottom: 5px;", "Plot 3"),
          div(id="plot3-container", class="plot-container", plotOutput("plot3", height="400px"))
        ),
        column(6, 
          div(style="text-align: center; font-weight: bold; color: #2c3e50; margin-bottom: 5px;", "Plot 4"),
          div(id="plot4-container", class="plot-container", plotOutput("plot4", height="400px"))
        )
      )
    )
  )),
  
  tabPanel("Differential Expression",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("DE Parameters"),
        selectInput("de_group_by", "Comparison Group", choices = NULL),
        fluidRow(
          column(6, selectInput("de_ident_1", "Ident 1", choices = NULL)),
          column(6, selectInput("de_ident_2", "Ident 2", choices = NULL))
        ),
        p(tags$small("Select 'All Others' in Ident 2 for one-vs-rest comparison.")),
        hr(),
        numericInput("de_logfc", "LogFC Threshold", value = 0.25, min = 0, max = 5, step = 0.05),
        numericInput("de_minpct", "Min Pct", value = 0.1, min = 0, max = 1, step = 0.05),
        numericInput("de_pval", "P-value Threshold", value = 0.05, min = 0.001, max = 0.1, step = 0.01),
        selectInput("de_test", "Test Type", choices = c("wilcox", "t-test", "roc", "LR"), selected = "wilcox"),
        br(),
        actionButton("run_de", "Run Differential Expression", class="btn-primary btn-block"),
        hr(),
        h5("Volcano Plot Settings"),
        sliderInput("volcano_pt_size", "Point Size", min = 0.5, max = 5, value = 2, step = 0.5),
        sliderInput("volcano_alpha", "Point Transparency", min = 0.1, max = 1, value = 0.6, step = 0.1),
        colourInput("volcano_sig_color", "Significant Color", value = "#E74C3C"),
        colourInput("volcano_nonsig_color", "Non-significant Color", value = "#95A5A6"),
        hr(),
        h5("Export"),
        selectInput("volcano_format", "Format", choices = c("png", "pdf", "jpg"), selected = "png"),
        numericInput("volcano_width", "Width (inches)", value = 10, min = 4, max = 20),
        numericInput("volcano_height", "Height (inches)", value = 8, min = 4, max = 20),
        downloadButton("download_volcano", "Download Volcano Plot", class="btn-info btn-block"),
        hr(),
        downloadButton("download_de", "Download CSV", class="btn-success btn-block")
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel("Results Table", 
            br(),
            DT::DTOutput("de_table")
          ),
          tabPanel("Volcano Plot", 
            br(),
            plotlyOutput("de_volcano", height = "600px")
          )
        )
      )
    )
  ),
  
  tabPanel("Pathway Enrichment",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Enrichment Analysis"),
        
        radioButtons("enrich_source", "Input Source",
                     choices = c("From DE Results" = "de", "Custom Gene List" = "custom"),
                     selected = "de"),
        
        conditionalPanel(
          condition = "input.enrich_source == 'custom'",
          textAreaInput("enrich_genes", "Gene List (one per line)", 
                       rows = 5, placeholder = "CD3D\nCD4\nCD8A\n...")
        ),
        
        hr(),
        h5("Analysis Settings"),
        
        selectInput("enrich_type", "Analysis Type",
                   choices = c("Over-Representation (ORA)" = "ora", 
                              "Gene Set Enrichment (GSEA)" = "gsea"),
                   selected = "ora"),
        
        conditionalPanel(
          condition = "input.enrich_type == 'ora' && input.enrich_source == 'de'",
          selectInput("enrich_gene_direction", "Gene Direction",
                     choices = c("All Significant Genes" = "all",
                                "Upregulated Only (positive log2FC)" = "up",
                                "Downregulated Only (negative log2FC)" = "down"),
                     selected = "all")
        ),
        
        selectInput("enrich_database", "Database",
                   choices = c("MSigDB Hallmarks" = "hallmark",
                              "GO Biological Process" = "go_bp",
                              "GO Molecular Function" = "go_mf",
                              "GO Cellular Component" = "go_cc",
                              "KEGG" = "kegg",
                              "Reactome" = "reactome"),
                   selected = "hallmark"),
        
        selectInput("enrich_organism", "Organism",
                   choices = c("Human" = "human", "Mouse" = "mouse"),
                   selected = "human"),
        
        hr(),
        h5("Parameters"),
        
        numericInput("enrich_pval", "P-value Cutoff", value = 0.05, min = 0.001, max = 0.1, step = 0.01),
        numericInput("enrich_qval", "Q-value Cutoff", value = 0.05, min = 0.001, max = 0.1, step = 0.01),
        numericInput("enrich_minsize", "Min Gene Set Size", value = 10, min = 5, max = 100),
        numericInput("enrich_maxsize", "Max Gene Set Size", value = 500, min = 100, max = 2000),
        
        br(),
        actionButton("run_enrichment", "Run Enrichment Analysis", class="btn-primary btn-block"),
        
        hr(),
        h5("Export Results"),
        downloadButton("download_enrichment", "Download Results (CSV)", class="btn-success btn-block"),
        
        hr(),
        h5("Export Plots"),
        selectInput("enrich_plot_format", "Format", choices = c("png", "pdf", "jpg"), selected = "png"),
        numericInput("enrich_plot_width", "Width (inches)", value = 10, min = 4, max = 20),
        numericInput("enrich_plot_height", "Height (inches)", value = 8, min = 4, max = 20),
        
        conditionalPanel(
          condition = "input.enrich_type == 'ora'",
          downloadButton("download_enrichment_dotplot", "Download Dot Plot", class="btn-info btn-block"),
          downloadButton("download_enrichment_barplot", "Download Bar Plot", class="btn-info btn-block"),
          downloadButton("download_enrichment_network", "Download Network Plot", class="btn-info btn-block")
        ),
        
        conditionalPanel(
          condition = "input.enrich_type == 'gsea'",
          downloadButton("download_enrichment_gsea", "Download GSEA Curves", class="btn-info btn-block")
        )
      ),
      mainPanel(
        width = 9,
        uiOutput("enrichment_plots_ui")
      )
    )
  ),
  
  tabPanel("Heatmap",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Heatmap Settings"),
        selectizeInput("hm_features", "Features (select multiple)", choices = NULL, multiple = TRUE, options = list(
          placeholder = 'Type genes or metadata...'
        )),
        selectInput("hm_group_by", "Group Cells By", choices = NULL),
        selectInput("hm_scaling", "Scaling / Z-score", 
                    choices = c("None (Real Expression)" = "none", 
                                "Row Z-score" = "row")),
        checkboxInput("hm_raster", "Rasterize", value = TRUE),
        hr(),
        h5("Annotation Colors"),
        selectInput("hm_anno_color_source", "Color Source", choices = c("Default", "Palette", "Manual")),
        conditionalPanel(
          condition = "input.hm_anno_color_source == 'Palette'",
          selectInput("hm_anno_palette", "Annotation Palette", 
                      choices = list(
                        "Viridis" = c("viridis", "magma", "plasma", "inferno", "cividis"),
                        "RColorBrewer" = c("Set1", "Set2", "Set3", "Dark2", "Paired", "Pastel1", "Accent"),
                        "MetBrewer" = names(MetBrewer::MetPalettes)
                      ))
        ),
        conditionalPanel(
          condition = "input.hm_anno_color_source == 'Manual'",
          uiOutput("hm_manual_color_ui")
        ),
        
        hr(),
        h5("Dimensions & Aesthetics"),
        fluidRow(
          column(6, numericInput("hm_width", "Width (in)", value = 12, min = 1, max = 50)),
          column(6, numericInput("hm_height", "Height (in)", value = 14, min = 1, max = 50))
        ),
        numericInput("hm_lines_width", "Gap Width", value = 1, min = 0, max = 10, step = 0.1),
        selectInput("hm_palette", "Color Palette", 
                    choices = c("Viridis" = "viridis", "Magma" = "magma", "Plasma" = "plasma", "Inferno" = "inferno", "Cividis" = "cividis", "Spectral" = "Spectral", "RdBu" = "RdBu", "RdYlBu" = "RdYlBu")),
        numericInput("hm_font_size", "Gene Font Size", value = 10, min = 1, max = 20),
        numericInput("hm_legend_font_size", "Legend Font Size", value = 10, min = 1, max = 20),
        br(),
        actionButton("run_hm", "Generate Heatmap", class = "btn-primary btn-block"),
        hr(),
        downloadButton("download_hm", "Download Heatmap", class = "btn-success btn-block")
      ),
      mainPanel(
        width = 9,
        plotOutput("heatmap_plot", height = "800px")
      )
    )
  )
)

server <- function(input, output, session) {
  
  original_obj <- reactiveVal(NULL)
  seurat_obj <- reactiveVal(NULL)
  
  observeEvent(input$seurat_file, {
    req(input$seurat_file)
    tryCatch({
      path <- input$seurat_file$datapath
      ext <- tolower(tools::file_ext(input$seurat_file$name))
      
      obj <- NULL
      if (ext == "rds") {
        obj <- readRDS(path)
      } else if (ext == "h5ad") {
        # Strategy 1: SeuratDisk (Fast, but fragile)
        if (requireNamespace("SeuratDisk", quietly = TRUE)) {
          showNotification("Attempting conversion with SeuratDisk...", type="message", duration=5)
          dest <- paste0(path, ".h5seurat")
          tryCatch({
            SeuratDisk::Convert(path, dest = dest, overwrite = TRUE)
            obj <- SeuratDisk::LoadH5Seurat(dest)
          }, error = function(e) {
            showNotification(paste("SeuratDisk failed:", e$message, "Trying fallback..."), type="warning")
            obj <<- NULL # Signal failure to fallback
          }, finally = {
            if(file.exists(dest)) unlink(dest)
          })
        }
        
        # Strategy 2: Zellkonverter (Robust, but requires Python env setup)
        if (is.null(obj)) {
          if (requireNamespace("zellkonverter", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
            showNotification("Attempting conversion with zellkonverter (this may take a while first time)...", type="message", duration=10)
            tryCatch({
              sce <- zellkonverter::readH5AD(path)
              obj <- as.Seurat(sce, counts = "X", data = "X")
            }, error = function(e) {
              stop(paste("All loading methods failed. Zellkonverter error:", e$message))
            })
          } else {
            stop("SeuratDisk failed and zellkonverter is not installed for fallback.")
          }
        }
      } else {
        stop("Unsupported file extension")
      }
      
      # --- Seurat Object Normalization ---
      # Fix Assay Names (zellkonverter often uses 'originalexp')
      if ("originalexp" %in% Assays(obj) && !"RNA" %in% Assays(obj)) {
        obj <- RenameAssays(obj, originalexp = "RNA")
      }
      # If 'RNA' is still not default (e.g. it was 'Spatial' or something else), force it if 'RNA' exists
      if ("RNA" %in% Assays(obj)) {
        DefaultAssay(obj) <- "RNA"
      }
      
      # Fix Reduction Names (zellkonverter often uses 'X_umap', 'X_pca')
      for (red in names(obj@reductions)) {
        if (red == "X_umap") {
          obj@reductions$umap <- obj@reductions$X_umap
          obj@reductions$X_umap <- NULL
          Key(obj@reductions$umap) <- "umap_"
        } else if (red == "X_pca") {
          obj@reductions$pca <- obj@reductions$X_pca
          obj@reductions$X_pca <- NULL
          Key(obj@reductions$pca) <- "PC_"
        } else if (red == "X_tsne") {
          obj@reductions$tsne <- obj@reductions$X_tsne
          obj@reductions$X_tsne <- NULL
          Key(obj@reductions$tsne) <- "tSNE_"
        }
      }
      # -----------------------------------
      
      original_obj(obj)
      seurat_obj(obj)
      showNotification("Loaded successfully!", type="message", duration=3)
    }, error = function(e) showNotification(paste("Error:", e$message), type="error"))
  })
  
  # --- Convert Ensembl IDs ---
  observeEvent(input$convert_ensembl, {
    req(seurat_obj())
    
    species <- input$ensembl_species
    db_pkg <- if(species == "Human") "org.Hs.eg.db" else "org.Mm.eg.db"
    
    if (!requireNamespace(db_pkg, quietly = TRUE)) {
      showNotification(paste(db_pkg, "is not installed."), type="error")
      return()
    }
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      showNotification("AnnotationDbi is not installed.", type="error")
      return()
    }
    
    # Load the package dynamically
    library(db_pkg, character.only = TRUE)
    # Get the db object (usually named same as package)
    org_db <- get(db_pkg)
    
    withProgress(message = "Converting IDs...", value = 0, {
      tryCatch({
        obj <- seurat_obj()
        # Get current feature names
        current_ids <- rownames(obj)
        
        # Map IDs
        incProgress(0.2, detail = paste("Mapping", species, "Ensembl to Symbols"))
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
          showNotification("No Ensembl IDs mapped to Symbols.", type="warning")
          return()
        }
        
        # Create new names vector: Use Symbol if mapped, else keep original ID
        new_names <- current_ids
        new_names[mapped_indices] <- mapped_symbols[mapped_indices]
        
        # Handle duplicates (e.g. distinct Ensembl IDs mapping to same Symbol)
        # We make them unique by appending .1, .2 etc.
        incProgress(0.5, detail = "Handling duplicates")
        new_names <- make.unique(new_names)
        
        # Rename features in the object
        # Since RenameGenesSeurat doesn't exist, we create a new assay with swapped counts
        incProgress(0.7, detail = "Updating Seurat Object")
        
        # Get counts matrix
        # Use LayerData if available (Seurat v5), or GetAssayData
        counts_mat <- GetAssayData(obj, layer = "counts")
        
        # Rename rows
        rownames(counts_mat) <- new_names
        
        # Create new assay
        new_assay <- CreateAssayObject(counts = counts_mat)
        
        # Try to bring over 'data' layer if it exists and has same structure
        if ("data" %in% Layers(obj, assay = DefaultAssay(obj))) {
           data_mat <- GetAssayData(obj, layer = "data")
           if (nrow(data_mat) == nrow(counts_mat)) {
             rownames(data_mat) <- new_names
             new_assay <- SetAssayData(new_assay, layer = "data", new.data = data_mat)
           }
        }
        
        # Replace the default assay (usually RNA) with our renamed version
        # It's safer to replace the assay entirely to ensure consistency
        assay_name <- DefaultAssay(obj)
        obj[[assay_name]] <- new_assay
        
        # If there were reductions, they might still have old feature names in loadings
        # For now we keep them but note that feature loadings might be mismatched.
        # Ideally we'd rename rows of Loadings(obj, reduction) as well.
        for (red in Reductions(obj)) {
           loadings <- Loadings(obj, reduction = red)
           if (!is.null(loadings) && nrow(loadings) == length(current_ids)) {
             # Only rename if dimensions match exactly (meaning all features were used)
             # Often PCA uses only variable features, so this might not match exactly.
             # Safe skip for now or we could try matching names.
           }
        }
        
        seurat_obj(obj)
        
        # Force UI updates
        # Ensure all feature dropdowns are updated with new symbols
        all_features <- c(rownames(obj), colnames(obj@meta.data))
        lapply(c("p1","p2","p3","p4"), function(id) {
          updateSelectizeInput(session, paste0(id, "-feature"), choices=all_features, server = TRUE)
        })
        
        showNotification(paste("Renamed", num_mapped, "features to Symbols!"), type="message")
        
      }, error = function(e) {
        showNotification(paste("Conversion Error:", e$message), type="error")
      })
    })
  })
  
  # --- DE Reactive Values & Logic ---
  de_results <- reactiveVal(NULL)
  
  # Update DE metadata column choices
  observe({
    req(seurat_obj())
    meta_cols <- names(seurat_obj()@meta.data)
    updateSelectInput(session, "de_group_by", choices = meta_cols, selected = "seurat_clusters")
  })
  
  # Update Ident 1 and Ident 2 based on chosen metadata column
  observeEvent(input$de_group_by, {
    req(seurat_obj(), input$de_group_by)
    lvls <- sort(unique(as.character(seurat_obj()@meta.data[[input$de_group_by]])))
    updateSelectInput(session, "de_ident_1", choices = lvls)
    updateSelectInput(session, "de_ident_2", choices = c("All Others", lvls))
  })
  
  observeEvent(input$run_de, {
    req(seurat_obj(), input$de_group_by, input$de_ident_1)
    
    withProgress(message = "Calculating Markers...", value = 0, {
      tryCatch({
        obj <- seurat_obj()
        # Set Idents safely on a copy or isolated
        Idents(obj) <- input$de_group_by
        
        ident2 <- if(input$de_ident_2 == "All Others") NULL else input$de_ident_2
        
        incProgress(0.5, detail = "Running FindMarkers (this may take a minute)")
        
        res <- FindMarkers(
          obj, 
          ident.1 = input$de_ident_1, 
          ident.2 = ident2,
          logfc.threshold = input$de_logfc,
          min.pct = input$de_minpct,
          test.use = input$de_test
        )
        
        res$gene <- rownames(res)
        de_results(res)
        
        showNotification("DE complete!", type="message")
      }, error = function(e) {
        showNotification(paste("DE Error:", e$message), type="error")
      })
    })
  })
  
  output$de_table <- DT::renderDT({
    req(de_results())
    DT::datatable(de_results(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) %>%
      DT::formatRound(columns = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), digits = 4)
  })
  
  # Reactive volcano plot (for display and download)
  volcano_plot <- reactive({
    req(de_results())
    df <- de_results()
    
    # Get threshold values from inputs
    logfc_thresh <- input$de_logfc
    pval_thresh <- input$de_pval
    if (is.null(pval_thresh)) pval_thresh <- 0.05
    
    # Get aesthetic controls
    pt_size <- input$volcano_pt_size
    if (is.null(pt_size)) pt_size <- 2
    alpha_val <- input$volcano_alpha
    if (is.null(alpha_val)) alpha_val <- 0.6
    sig_color <- input$volcano_sig_color
    if (is.null(sig_color)) sig_color <- "#E74C3C"
    nonsig_color <- input$volcano_nonsig_color
    if (is.null(nonsig_color)) nonsig_color <- "#95A5A6"
    
    # Determine significance based on user's threshold
    df$significant <- df$p_val_adj < pval_thresh
    
    # Calculate symmetric x-axis limits
    max_abs_fc <- max(abs(df$avg_log2FC), na.rm = TRUE)
    x_limit <- ceiling(max_abs_fc * 1.1)  # Add 10% padding
    
    ggplot(df, aes(x = avg_log2FC, y = -log10(p_val), text = gene, color = significant)) +
      geom_point(size = pt_size, alpha = alpha_val) +
      # Add threshold lines
      geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "black", linewidth = 0.8) +
      geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "black", linewidth = 0.8) +
      # Symmetric x-axis
      coord_cartesian(xlim = c(-x_limit, x_limit)) +
      theme_minimal() +
      scale_color_manual(values = c(nonsig_color, sig_color), name = paste0("p-adj < ", pval_thresh)) +
      labs(title = paste("Volcano Plot:", input$de_ident_1, "vs", input$de_ident_2),
           x = "avg_log2FC", y = "-log10(p-value)")
  })
  
  output$de_volcano <- renderPlotly({
    req(volcano_plot())
    ggplotly(volcano_plot(), tooltip = "text")
  })
  
  output$download_volcano <- downloadHandler(
    filename = function() { 
      paste0("Volcano_", input$de_ident_1, "_vs_", input$de_ident_2, "_", Sys.Date(), ".", input$volcano_format) 
    },
    content = function(file) {
      req(volcano_plot())
      ggsave(file, volcano_plot(), 
             width = input$volcano_width, 
             height = input$volcano_height, 
             device = input$volcano_format)
    }
  )
  
  output$download_de <- downloadHandler(
    filename = function() { paste0("DE_", input$de_ident_1, "_vs_", input$de_ident_2, "_", Sys.Date(), ".csv") },
    content = function(file) {
      req(de_results())
      write.csv(de_results(), file, row.names = FALSE)
    }
  )
  # ===== PATHWAY ENRICHMENT LOGIC =====
  
  enrichment_results <- reactiveVal(NULL)
  
  # Helper function: Convert gene symbols to Entrez IDs
  convert_to_entrez <- function(genes, organism = "human") {
    if (!has_enrichment) return(NULL)
    
    tryCatch({
      if (organism == "human") {
        suppressMessages(library(org.Hs.eg.db))
        org_db <- org.Hs.eg.db
      } else {
        suppressMessages(library(org.Mm.eg.db))
        org_db <- org.Mm.eg.db
      }
      
      gene_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", 
                      OrgDb = org_db, drop = TRUE)
      return(gene_map)
    }, error = function(e) {
      showNotification(paste("Gene ID conversion error:", e$message), type = "error")
      return(NULL)
    })
  }
  
  # Helper function: Get gene sets from database
  get_gene_sets <- function(database, organism) {
    if (!has_enrichment) return(NULL)
    
    tryCatch({
      species <- ifelse(organism == "human", "Homo sapiens", "Mus musculus")
      
      if (database == "hallmark") {
        msigdbr(species = species, category = "H")
      } else if (database == "reactome") {
        msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME")
      } else {
        # For GO and KEGG, we'll use clusterProfiler's built-in functions
        NULL
      }
    }, error = function(e) {
      showNotification(paste("Database loading error:", e$message), type = "error")
      return(NULL)
    })
  }
  
  # Run enrichment analysis
  observeEvent(input$run_enrichment, {
    req(has_enrichment)
    
    # Get gene list
    if (input$enrich_source == "de") {
      req(de_results())
      df <- de_results()
      
      # For ORA: use only significant genes (typically upregulated/significant)
      # For GSEA: use all genes ranked by fold change
      if (input$enrich_type == "ora") {
        # Filter to significant genes (p_val_adj < 0.05)
        df_sig <- df[df$p_val_adj < 0.05, ]
        if (nrow(df_sig) == 0) {
          showNotification("No significant genes found for ORA. Please adjust DE p-value threshold.", 
                          type = "error")
          return(NULL)
        }
        
        # Apply direction filter
        direction <- input$enrich_gene_direction
        if (is.null(direction)) direction <- "all"
        
        if (direction == "up") {
          # Upregulated in Group 1 (positive log2FC)
          df_sig <- df_sig[df_sig$avg_log2FC > 0, ]
          direction_label <- "upregulated"
        } else if (direction == "down") {
          # Downregulated in Group 1 (negative log2FC)
          df_sig <- df_sig[df_sig$avg_log2FC < 0, ]
          direction_label <- "downregulated"
        } else {
          direction_label <- "significant"
        }
        
        if (nrow(df_sig) == 0) {
          showNotification(paste("No", direction_label, "genes found. Try different direction."), 
                          type = "error")
          return(NULL)
        }
        
        genes <- df_sig$gene
        gene_fc <- setNames(df_sig$avg_log2FC, df_sig$gene)
        showNotification(paste("Using", length(genes), direction_label, "genes from DE results"), 
                        type = "message")
      } else {
        # GSEA uses all genes ranked by fold change
        genes <- df$gene
        gene_fc <- setNames(df$avg_log2FC, df$gene)
        showNotification(paste("Using", length(genes), "ranked genes for GSEA"), 
                        type = "message")
      }
    } else {
      req(input$enrich_genes)
      genes <- strsplit(input$enrich_genes, "\\n")[[1]]
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      gene_fc <- NULL
    }
    
    if (length(genes) == 0) {
      showNotification("No genes provided", type = "error")
      return(NULL)
    }
    
    # Show progress
    withProgress(message = "Running enrichment analysis...", value = 0, {
      
      # Convert to Entrez IDs
      incProgress(0.2, detail = "Converting gene IDs")
      gene_map <- convert_to_entrez(genes, input$enrich_organism)
      
      if (is.null(gene_map) || nrow(gene_map) == 0) {
        showNotification("No genes could be mapped to Entrez IDs", type = "error")
        return(NULL)
      }
      
      showNotification(paste(nrow(gene_map), "out of", length(genes), "genes mapped"), 
                      type = "message")
      
      # Get organism database
      if (input$enrich_organism == "human") {
        suppressMessages(library(org.Hs.eg.db))
        org_db <- org.Hs.eg.db
      } else {
        suppressMessages(library(org.Mm.eg.db))
        org_db <- org.Mm.eg.db
      }
      
      # Run analysis based on type
      incProgress(0.3, detail = "Running enrichment")
      
      result <- tryCatch({
        if (input$enrich_type == "ora") {
          # Over-Representation Analysis
          if (input$enrich_database %in% c("go_bp", "go_mf", "go_cc")) {
            ont <- switch(input$enrich_database,
                         "go_bp" = "BP",
                         "go_mf" = "MF",
                         "go_cc" = "CC")
            enrichGO(gene = gene_map$ENTREZID,
                    OrgDb = org_db,
                    ont = ont,
                    pAdjustMethod = "BH",
                    pvalueCutoff = input$enrich_pval,
                    qvalueCutoff = input$enrich_qval,
                    minGSSize = input$enrich_minsize,
                    maxGSSize = input$enrich_maxsize)
          } else if (input$enrich_database == "kegg") {
            organism_code <- ifelse(input$enrich_organism == "human", "hsa", "mmu")
            enrichKEGG(gene = gene_map$ENTREZID,
                      organism = organism_code,
                      pAdjustMethod = "BH",
                      pvalueCutoff = input$enrich_pval,
                      qvalueCutoff = input$enrich_qval,
                      minGSSize = input$enrich_minsize,
                      maxGSSize = input$enrich_maxsize)
          } else {
            # MSigDB (Hallmark, Reactome)
            gene_sets <- get_gene_sets(input$enrich_database, input$enrich_organism)
            if (is.null(gene_sets)) {
              showNotification("Failed to load gene sets", type = "error")
              return(NULL)
            }
            
            term2gene <- gene_sets %>% dplyr::select(gs_name, entrez_gene)
            enricher(gene = gene_map$ENTREZID,
                    TERM2GENE = term2gene,
                    pAdjustMethod = "BH",
                    pvalueCutoff = input$enrich_pval,
                    qvalueCutoff = input$enrich_qval,
                    minGSSize = input$enrich_minsize,
                    maxGSSize = input$enrich_maxsize)
          }
        } else {
          # GSEA
          # Prepare ranked gene list
          if (is.null(gene_fc)) {
            showNotification("GSEA requires fold change values. Please use DE results.", 
                           type = "error")
            return(NULL)
          }
          
          # Map fold changes to Entrez IDs
          gene_list <- gene_fc[gene_map$SYMBOL]
          names(gene_list) <- gene_map$ENTREZID[match(names(gene_list), gene_map$SYMBOL)]
          gene_list <- sort(gene_list, decreasing = TRUE)
          gene_list <- gene_list[!is.na(names(gene_list))]
          
          if (input$enrich_database %in% c("go_bp", "go_mf", "go_cc")) {
            ont <- switch(input$enrich_database,
                         "go_bp" = "BP",
                         "go_mf" = "MF",
                         "go_cc" = "CC")
            gseGO(geneList = gene_list,
                 OrgDb = org_db,
                 ont = ont,
                 pAdjustMethod = "BH",
                 pvalueCutoff = input$enrich_pval,
                 minGSSize = input$enrich_minsize,
                 maxGSSize = input$enrich_maxsize)
          } else if (input$enrich_database == "kegg") {
            organism_code <- ifelse(input$enrich_organism == "human", "hsa", "mmu")
            gseKEGG(geneList = gene_list,
                   organism = organism_code,
                   pAdjustMethod = "BH",
                   pvalueCutoff = input$enrich_pval,
                   minGSSize = input$enrich_minsize,
                   maxGSSize = input$enrich_maxsize)
          } else {
            # MSigDB (Hallmark, Reactome)
            gene_sets <- get_gene_sets(input$enrich_database, input$enrich_organism)
            if (is.null(gene_sets)) {
              showNotification("Failed to load gene sets", type = "error")
              return(NULL)
            }
            
            term2gene <- gene_sets %>% dplyr::select(gs_name, entrez_gene)
            GSEA(geneList = gene_list,
                TERM2GENE = term2gene,
                pAdjustMethod = "BH",
                pvalueCutoff = input$enrich_pval,
                minGSSize = input$enrich_minsize,
                maxGSSize = input$enrich_maxsize)
          }
        }
      }, error = function(e) {
        showNotification(paste("Enrichment error:", e$message), type = "error")
        return(NULL)
      })
      
      incProgress(0.5, detail = "Processing results")
      
      if (!is.null(result) && nrow(result@result) > 0) {
        # Convert Entrez IDs to gene symbols for better readability
        if (input$enrich_organism == "human") {
          suppressMessages(library(org.Hs.eg.db))
          org_db <- org.Hs.eg.db
        } else {
          suppressMessages(library(org.Mm.eg.db))
          org_db <- org.Mm.eg.db
        }
        
        result <- setReadable(result, OrgDb = org_db, keyType = "ENTREZID")
        enrichment_results(result)
        showNotification(paste("Found", nrow(result@result), "enriched pathways"), 
                        type = "message")
      } else {
        showNotification("No significant enrichment found", type = "warning")
        enrichment_results(NULL)
      }
    })
  })
  
  # Dynamic UI for enrichment plots based on analysis type
  output$enrichment_plots_ui <- renderUI({
    if (input$enrich_type == "ora") {
      tabsetPanel(
        tabPanel("Results Table",
          br(),
          DT::DTOutput("enrichment_table")
        ),
        tabPanel("Dot Plot",
          br(),
          plotOutput("enrichment_dotplot", height = "600px")
        ),
        tabPanel("Bar Plot",
          br(),
          plotOutput("enrichment_barplot", height = "600px")
        ),
        tabPanel("Network Plot",
          br(),
          plotOutput("enrichment_network", height = "700px")
        )
      )
    } else {
      # GSEA - only show results table and enrichment score plots
      tabsetPanel(
        tabPanel("Results Table",
          br(),
          DT::DTOutput("enrichment_table")
        ),
        tabPanel("GSEA Enrichment Curves",
          br(),
          plotOutput("enrichment_gsea", height = "800px")
        )
      )
    }
  })
  
  
  
  # Results table
  output$enrichment_table <- DT::renderDT({
    req(enrichment_results())
    result_df <- as.data.frame(enrichment_results())
    
    DT::datatable(result_df, 
                 options = list(pageLength = 25, scrollX = TRUE),
                 filter = "top") %>%
      DT::formatRound(columns = c("pvalue", "p.adjust", "qvalue"), digits = 4)
  })
  
  # Reactive plot objects for reuse in display and download
  enrichment_dotplot_obj <- reactive({
    req(enrichment_results())
    dotplot(enrichment_results(), showCategory = 20) + 
      theme(axis.text.y = element_text(size = 10))
  })
  
  enrichment_barplot_obj <- reactive({
    req(enrichment_results())
    barplot(enrichment_results(), showCategory = 20) +
      theme(axis.text.y = element_text(size = 10))
  })
  
  enrichment_network_obj <- reactive({
    req(enrichment_results())
    result_sim <- pairwise_termsim(enrichment_results())
    emapplot(result_sim, showCategory = 30)
  })
  
  # Dot plot
  output$enrichment_dotplot <- renderPlot({
    req(enrichment_dotplot_obj())
    tryCatch({
      enrichment_dotplot_obj()
    }, error = function(e) {
      ggplot() + annotate("text", x = 0, y = 0, 
                         label = paste("Plot error:", e$message), 
                         size = 4, color = "red") + theme_void()
    })
  })
  
  
  # Bar plot
  output$enrichment_barplot <- renderPlot({
    req(enrichment_barplot_obj())
    tryCatch({
      enrichment_barplot_obj()
    }, error = function(e) {
      ggplot() + annotate("text", x = 0, y = 0, 
                         label = paste("Plot error:", e$message), 
                         size = 4, color = "red") + theme_void()
    })
  })
  
  
  # Network plot
  output$enrichment_network <- renderPlot({
    req(enrichment_network_obj())
    tryCatch({
      enrichment_network_obj()
    }, error = function(e) {
      ggplot() + annotate("text", x = 0, y = 0, 
                         label = paste("Network plot error:", e$message), 
                         size = 4, color = "red") + theme_void()
    })
  })
  
  # Reactive GSEA plot object for reuse
  enrichment_gsea_obj <- reactive({
    req(enrichment_results(), input$enrich_type == "gsea")
    n_pathways <- min(5, nrow(enrichment_results()))
    gseaplot2(enrichment_results(), geneSetID = 1:n_pathways, 
             pvalue_table = TRUE)
  })
  
  # GSEA enrichment plots
  output$enrichment_gsea <- renderPlot({
    req(enrichment_gsea_obj())
    tryCatch({
      enrichment_gsea_obj()
    }, error = function(e) {
      ggplot() + annotate("text", x = 0, y = 0, 
                         label = paste("GSEA plot error:", e$message), 
                         size = 4, color = "red") + theme_void()
    })
  })
  
  
  # Download enrichment results
  output$download_enrichment <- downloadHandler(
    filename = function() { 
      paste0("Enrichment_", input$enrich_database, "_", Sys.Date(), ".csv") 
    },
    content = function(file) {
      req(enrichment_results())
      write.csv(as.data.frame(enrichment_results()), file, row.names = FALSE)
    }
  )
  
  # Download dot plot
  output$download_enrichment_dotplot <- downloadHandler(
    filename = function() { 
      paste0("Enrichment_Dotplot_", input$enrich_database, "_", Sys.Date(), ".", input$enrich_plot_format) 
    },
    content = function(file) {
      req(enrichment_dotplot_obj())
      ggsave(file, enrichment_dotplot_obj(), 
             width = input$enrich_plot_width, 
             height = input$enrich_plot_height, 
             device = input$enrich_plot_format)
    }
  )
  
  # Download bar plot
  output$download_enrichment_barplot <- downloadHandler(
    filename = function() { 
      paste0("Enrichment_Barplot_", input$enrich_database, "_", Sys.Date(), ".", input$enrich_plot_format) 
    },
    content = function(file) {
      req(enrichment_barplot_obj())
      ggsave(file, enrichment_barplot_obj(), 
             width = input$enrich_plot_width, 
             height = input$enrich_plot_height, 
             device = input$enrich_plot_format)
    }
  )
  
  # Download network plot
  output$download_enrichment_network <- downloadHandler(
    filename = function() { 
      paste0("Enrichment_Network_", input$enrich_database, "_", Sys.Date(), ".", input$enrich_plot_format) 
    },
    content = function(file) {
      req(enrichment_network_obj())
      ggsave(file, enrichment_network_obj(), 
             width = input$enrich_plot_width, 
             height = input$enrich_plot_height, 
             device = input$enrich_plot_format)
    }
  )
  
  # Download GSEA enrichment curves
  output$download_enrichment_gsea <- downloadHandler(
    filename = function() { 
      paste0("Enrichment_GSEA_", input$enrich_database, "_", Sys.Date(), ".", input$enrich_plot_format) 
    },
    content = function(file) {
      req(enrichment_gsea_obj())
      ggsave(file, enrichment_gsea_obj(), 
             width = input$enrich_plot_width, 
             height = input$enrich_plot_height, 
             device = input$enrich_plot_format)
    }
  )
  
  # ===== END PATHWAY ENRICHMENT LOGIC =====
  
  output$obj_info <- renderText({ 
    req(seurat_obj())
    paste(ncol(seurat_obj()), "cells,", nrow(seurat_obj()), "features")
  })
  
  # Preview plot in subsetting tab
  output$preview_plot <- renderPlot({
    req(seurat_obj())
    tryCatch({
      red <- names(seurat_obj()@reductions)[1]
      if (!is.null(red)) {
        DimPlot(seurat_obj(), reduction=red, pt.size=0.5) + 
          theme(legend.position="right", legend.text=element_text(size=8))
      } else {
        ggplot() + annotate("text", x=0, y=0, label="No reductions available") + theme_void()
      }
    }, error=function(e) ggplot() + theme_void())
  })
  
  # --- Subsetting ---
  observe({ 
    req(seurat_obj())
    updateSelectInput(session, "subset_col", choices=names(seurat_obj()@meta.data))
  })
  
  observeEvent(input$subset_col, { 
    req(seurat_obj())
    vals <- sort(unique(seurat_obj()@meta.data[[input$subset_col]]))
    updateSelectizeInput(session, "subset_levels", choices=vals)
  })
  
  observeEvent(input$apply_subset_meta, {
    req(seurat_obj(), input$subset_levels)
    cells <- rownames(seurat_obj()@meta.data[seurat_obj()@meta.data[[input$subset_col]] %in% input$subset_levels, ])
    seurat_obj(subset(seurat_obj(), cells=cells))
    showNotification(paste("Filtered to", length(cells), "cells"), type="message")
  })
  
  observeEvent(input$reset_data, { 
    seurat_obj(original_obj())
    showNotification("Reset to original", type="message")
  })
  
  # --- UCell ---
  observeEvent(input$calc_ucell, {
    req(seurat_obj(), input$sig_genes, input$sig_name)
    if(!has_ucell) { 
      showNotification("UCell not installed", type="error")
      return()
    }
    
    genes <- trimws(unlist(strsplit(input$sig_genes, "\n")))
    genes <- genes[genes != ""]
    
    if(length(genes) == 0) {
      showNotification("No genes provided", type="warning")
      return()
    }
    
    sig_list <- list()
    sig_list[[input$sig_name]] <- genes
    
    tryCatch({
      # Use isolate to prevent triggering plot updates
      obj <- isolate(seurat_obj())
      obj <- AddModuleScore_UCell(obj, features=sig_list, name=NULL)
      
      # Find the actual column name
      sig_col <- grep(input$sig_name, colnames(obj@meta.data), value=TRUE)[1]
      
      # Update object without triggering plots
      isolate({
        seurat_obj(obj)
        
        # Force update of all feature dropdowns
        all_features <- c(rownames(obj), colnames(obj@meta.data))
        lapply(c("p1","p2","p3","p4"), function(id) {
          updateSelectizeInput(session, paste0(id, "-feature"), choices=all_features)
        })
      })
      
      showNotification(paste("UCell score added:", sig_col, "\nYou can now plot it in FeaturePlot or ViolinPlot"), type="message", duration=8)
      
    }, error=function(e) {
      showNotification(paste("UCell error:", e$message), type="error")
    })
  })
  
  # --- Color Helper (Refactored to use color_utils.R) ---
  get_colors <- function(id, obj, group_var, for_scpubr = FALSE) {
    ns <- function(x) paste0(id, "-", x)
    
    # Get levels for the grouping variable
    levels <- get_group_levels(obj, group_var)
    
    # Use unified color utility function
    colors <- get_plot_colors(input, ns, levels, named = for_scpubr)
    
    return(colors)
  }
  
  # --- Plot Generator ---
  generate_plot <- function(id, obj) {
    ns <- function(x) paste0(id, "-", x)
    ptype <- input[[ns("plot_type")]]
    req(ptype)
    
    red <- input[[ns("reduction")]]
    feat <- input[[ns("feature")]]
    grp <- input[[ns("group_by")]]
    splt <- input[[ns("split_by")]]
    
    if (!is.null(grp) && grp == "Default") grp <- NULL
    if (!is.null(splt) && splt == "None") splt <- NULL
    
    # Determine the correct grouping variable for color retrieval
    # For ClusterDistrBar, use cdb_group2 (fill variable) instead of group_by
    color_group_var <- grp
    if (ptype == "ClusterDistrBar") {
      color_group_var <- input[[ns("cdb_group2")]]
    }
    
    colors <- get_colors(id, obj, if(is.null(color_group_var)) "Default" else color_group_var)
    # For SCpubr, get named colors
    style <- input[[ns("plot_style")]]
    if (!is.null(style) && style == "SCpubr") {
      colors <- get_colors(id, obj, if(is.null(color_group_var)) "Default" else color_group_var, for_scpubr = TRUE)
    }
    pt_size <- input[[ns("pt_size")]]
    if (is.null(pt_size)) pt_size <- 1  # Default if not set
    
    p <- NULL
    
    if (ptype == "ClusterDistrBar") {
      # Use modular plotting function
      params <- validate_cluster_distribution_params(input, ns, obj)
      
      if (!is.null(params)) {
        p <- plot_cluster_distribution(
          obj = obj,
          group1 = params$group1,
          group2 = params$group2,
          colors = colors,
          flip = isTRUE(params$flip),
          show_counts = isTRUE(params$show_counts),
          count_size = if(is.null(params$count_size)) 3 else params$count_size,
          title = params$title,
          title_size = if(is.null(params$title_size)) 16 else params$title_size,
          show_legend = isTRUE(params$show_legend)
        )
      }
      
    } else {
      # Continuous colors for FeaturePlot
      cols_cont <- if(!is.null(colors) && length(colors) >= 2) c(colors[1], colors[length(colors)]) else NULL
      
      # Get plot style
      style <- input[[ns("plot_style")]]
      if (is.null(style)) style <- "Standard"
      
      # SCpubr plotting
      if (style == "SCpubr" && has_scpubr) {
        # Get SCpubr settings
        font_size <- input[[ns("scpubr_font_size")]]
        if (is.null(font_size)) font_size <- 14
        border <- input[[ns("scpubr_border")]]
        if (is.null(border)) border <- TRUE
        legend_pos <- input[[ns("scpubr_legend_pos")]]
        if (is.null(legend_pos)) legend_pos <- "bottom"  # SCpubr default
        raster_dpi <- input[[ns("scpubr_raster_dpi")]]
        if (is.null(raster_dpi)) raster_dpi <- 2048  # SCpubr default
        
        tryCatch({
          if (ptype == "DimPlot") {
            req(red)
            p <- SCpubr::do_DimPlot(
              sample = obj,
              reduction = red,
              group.by = grp,
              split.by = splt,
              colors.use = colors,
              pt.size = pt_size,
              font.size = font_size,
              plot_cell_borders = border,
              legend.position = legend_pos,
              raster = ncol(obj) > 50000,
              raster.dpi = raster_dpi
            )
          } else if (ptype == "FeaturePlot") {
            req(feat, red)
            # FeaturePlot - support viridis palettes
            use_viridis <- FALSE
            viridis_pal <- "D"  # Default to standard viridis
            pal_source <- input[[ns("color_source")]]
            if (!is.null(pal_source) && pal_source == "Palette") {
              pal_name <- input[[ns("palette_name")]]
              if (!is.null(pal_name) && pal_name %in% c("viridis", "magma", "plasma", "inferno", "cividis")) {
                use_viridis <- TRUE
                # Map to SCpubr's viridis.palette options (A-E)
                viridis_map <- list(magma="A", inferno="B", plasma="C", viridis="D", cividis="E")
                viridis_pal <- viridis_map[[pal_name]]
              }
            }
            
            p <- SCpubr::do_FeaturePlot(
              sample = obj,
              features = feat,
              reduction = red,
              split.by = splt,
              pt.size = pt_size,
              font.size = font_size,
              plot_cell_borders = border,
              legend.position = legend_pos,
              raster = ncol(obj) > 50000,
              raster.dpi = raster_dpi,
              use_viridis = use_viridis,
              viridis.palette = viridis_pal,
              viridis.direction = as.numeric(input[[ns("scpubr_viridis_dir")]])
            )
          } else if (ptype == "ViolinPlot") {
            req(feat)
            p <- tryCatch({
              result <- SCpubr::do_ViolinPlot(
                sample = obj,
                features = feat,
                group.by = grp,
                split.by = splt,
                colors.use = colors,
                font.size = font_size,
                legend.position = legend_pos
              )
            }, error = function(e) {
              ggplot() + annotate("text", x=0, y=0, label=paste("SCpubr Violin Error:", e$message), size=4, color="red") + theme_void()
            })
          } else if (ptype == "DotPlot") {
            req(feat)
            # DotPlot - support viridis palettes
            use_viridis <- FALSE
            viridis_pal <- "E"  # Default to cividis
            pal_source <- input[[ns("color_source")]]
            if (!is.null(pal_source) && pal_source == "Palette") {
              pal_name <- input[[ns("palette_name")]]
              if (!is.null(pal_name) && pal_name %in% c("viridis", "magma", "plasma", "inferno", "cividis")) {
                use_viridis <- TRUE
                # Map to SCpubr's viridis.palette options (A-E)
                viridis_map <- list(magma="A", inferno="B", plasma="C", viridis="D", cividis="E")
                viridis_pal <- viridis_map[[pal_name]]
              }
            }
            
            p <- tryCatch({
              result <- SCpubr::do_DotPlot(
                sample = obj,
                features = feat,
                group.by = grp,
                split.by = splt,
                font.size = font_size,
                legend.position = legend_pos,
                use_viridis = use_viridis,
                viridis.palette = viridis_pal,
                viridis.direction = as.numeric(input[[ns("scpubr_viridis_dir")]])
              )
            }, error = function(e) {
              # If split.by fails, try without it
              result <- SCpubr::do_DotPlot(
                sample = obj,
                features = feat,
                group.by = grp,
                font.size = font_size,
                legend.position = legend_pos,
                use_viridis = use_viridis,
                viridis.palette = viridis_pal,
                viridis.direction = as.numeric(input[[ns("scpubr_viridis_dir")]])
              )
            })
          }
          
          # Fix SCpubr Legend Orientation
          if (!is.null(p)) {
            leg_dir <- if (legend_pos %in% c("bottom", "top")) "horizontal" else "vertical"
            
            if (ptype == "FeaturePlot") {
              p <- p + guides(
                color = guide_colorbar(direction = leg_dir, title.position = "top", order = 1),
                fill = guide_colorbar(direction = leg_dir, title.position = "top", order = 1)
              )
            } else if (ptype == "ViolinPlot") {
              p <- p + guides(fill = guide_legend(direction = leg_dir, override.aes = list(size = 4), title.position = "top"))
            } else if (ptype == "DimPlot") {
              p <- p + guides(color = guide_legend(direction = leg_dir, override.aes = list(size = 4), title.position = "top"))
            } else if (ptype == "DotPlot") {
              p <- p + guides(
                color = guide_colorbar(direction = leg_dir, title.position = "top", order = 1),
                size = guide_legend(direction = leg_dir, title.position = "top", order = 2)
              )
            }
          }
          
        }, error = function(e) {
          # Fallback to standard if SCpubr fails
          p <<- ggplot() + annotate("text", x=0, y=0, label=paste("SCpubr error:", e$message), size=4, color="red") + theme_void()
        })
        
      } else {
        # Standard Seurat plotting
        if (ptype == "DimPlot") {
          req(red)  # Must have reduction
          p <- DimPlot(obj, reduction=red, group.by=grp, split.by=splt, cols=colors, pt.size=pt_size)
        } else if (ptype == "FeaturePlot") {
          req(feat, red)  # Must have feature AND reduction
          # Only pass cols if not NULL - let Seurat use default gradient otherwise
          if (!is.null(cols_cont)) {
            p <- FeaturePlot(obj, features=feat, reduction=red, split.by=splt, pt.size=pt_size, cols=cols_cont)
          } else {
            p <- FeaturePlot(obj, features=feat, reduction=red, split.by=splt, pt.size=pt_size)
          }
        } else if (ptype == "ViolinPlot") {
          req(feat)  # Must have feature
          p <- VlnPlot(obj, features=feat, group.by=grp, split.by=splt, cols=colors, pt.size=pt_size)
        } else if (ptype == "DotPlot") {
          req(feat)  # Must have feature(s)
          feats <- feat
          
          if (!is.null(cols_cont)) {
            p <- DotPlot(obj, features=feats, group.by=grp, split.by=splt, cols=cols_cont)
          } else {
            p <- DotPlot(obj, features=feats, group.by=grp, split.by=splt)
          }
        }
      }
    }
    
    # If p is still NULL, return a message
    if (is.null(p)) {
      return(ggplot() + annotate("text", x=0, y=0, label="Please configure plot settings", size=5, color="gray50") + theme_void())
    }
    
    # Customization
    title <- input[[ns("custom_title")]]
    if(isTruthy(title)) p <- p + ggtitle(title)
    
    # Flip coordinates if requested
    flip <- input[[ns("flip_coords")]]
    if (!is.null(flip) && flip && ptype %in% c("DotPlot", "ViolinPlot", "ClusterDistrBar")) {
      p <- p + coord_flip()
    }
    
    p <- p + theme(
      plot.title = element_text(size=input[[ns("title_size")]], hjust=as.numeric(input[[ns("title_align")]])),
      legend.position = if(input[[ns("show_legend")]]) "right" else "none"
    )
    
    return(p)
  }
  
  # --- Dynamic UI ---
  lapply(c("p1","p2","p3","p4"), function(id) {
    ns <- function(x) paste0(id, "-", x)
    
    output[[ns("dynamic_ui")]] <- renderUI({
      req(seurat_obj())
      ptype <- input[[ns("plot_type")]]
      reds <- names(seurat_obj()@reductions)
      metas <- c("Default", names(seurat_obj()@meta.data))
      
      # Get all features (genes + metadata for UCell)
      all_features <- c(rownames(seurat_obj()), colnames(seurat_obj()@meta.data))
      
      if (ptype == "ClusterDistrBar") {
        tagList(
          selectInput(ns("cdb_group1"), "X Axis (Sample)", choices=metas[-1]),
          selectInput(ns("cdb_group2"), "Fill (Cluster)", choices=metas[-1])
        )
      } else {
        tagList(
          if(ptype %in% c("DimPlot","FeaturePlot")) selectInput(ns("reduction"), "Reduction", choices=reds),
          if(ptype %in% c("DimPlot","ViolinPlot","DotPlot")) selectInput(ns("group_by"), "Group By", choices=metas),
          selectInput(ns("split_by"), "Split By", choices=c("None", metas[-1])),
          if(ptype == "DotPlot")
            selectizeInput(ns("feature"), "Features (select multiple)", choices=NULL, multiple=TRUE, options=list(
              placeholder = 'Type to search and select multiple genes...',
              onInitialize = I('function() { this.setValue(""); }')
            ))
          else if(ptype %in% c("FeaturePlot","ViolinPlot")) 
            selectizeInput(ns("feature"), "Feature", choices=NULL, options=list(
              placeholder = 'Type to search genes or metadata...',
              onInitialize = I('function() { this.setValue(""); }')
            ))
        )
      }
    })
    
    # Update feature choices server-side for autocomplete
    observe({
      req(seurat_obj())
      ptype <- input[[ns("plot_type")]]
      if (ptype %in% c("FeaturePlot","ViolinPlot","DotPlot")) {
        all_features <- c(rownames(seurat_obj()), colnames(seurat_obj()@meta.data))
        updateSelectizeInput(session, ns("feature"), choices=all_features, server = TRUE)
      }
    })
    
    # Point size UI - hide for ClusterDistrBar
    output[[ns("pt_size_ui")]] <- renderUI({
      ptype <- input[[ns("plot_type")]]
      if (is.null(ptype) || ptype != "ClusterDistrBar") {
        sliderInput(ns("pt_size"), "Point Size", min=0.1, max=5, value=1, step=0.1)
      }
    })
    
    # Orientation UI - show for DotPlot, ViolinPlot, ClusterDistrBar
    output[[ns("orientation_ui")]] <- renderUI({
      ptype <- input[[ns("plot_type")]]
      if (!is.null(ptype) && ptype %in% c("DotPlot", "ViolinPlot", "ClusterDistrBar")) {
        tagList(
          checkboxInput(ns("flip_coords"), "Flip to Horizontal", value=FALSE),
          if (ptype == "ClusterDistrBar") tagList(
            checkboxInput(ns("show_counts"), "Show Cell Counts", value=FALSE),
            conditionalPanel(
              condition = "input.show_counts", ns = ns,
              sliderInput(ns("count_size"), "Count Label Size", min=2, max=8, value=3, step=0.5)
            )
          )
        )
      }
    })
    
    # Plot style UI - hide SCpubr for ClusterDistrBar
    output[[ns("plot_style_ui")]] <- renderUI({
      ptype <- input[[ns("plot_type")]]
      if (is.null(ptype) || ptype == "ClusterDistrBar") {
        # No SCpubr for ClusterDistrBar
        selectInput(ns("plot_style"), "Plot Style", choices = c("Standard"))
      } else {
        selectInput(ns("plot_style"), "Plot Style", 
                    choices = c("Standard", if(has_scpubr) "SCpubr" else NULL))
      }
    })
    
    # Manual color UI
    output[[ns("manual_color_ui")]] <- renderUI({
      req(seurat_obj(), input[[ns("color_source")]] == "Manual")
      
      ptype <- input[[ns("plot_type")]]
      grp <- input[[ns("group_by")]]
      
      if (ptype == "ClusterDistrBar") grp <- input[[ns("cdb_group2")]]
      
      if (is.null(grp) || grp %in% c("Default", "None")) {
        lvls <- levels(Idents(seurat_obj()))
      } else {
        lvls <- sort(unique(seurat_obj()@meta.data[[grp]]))
      }
      
      # Create scrollable div for many levels
      color_pickers <- lapply(lvls, function(l) {
        colourInput(ns(paste0("col_", l)), paste("Color", l), value="gray")
      })
      
      if (length(lvls) > 20) {
        div(style="max-height: 400px; overflow-y: auto; padding: 5px;",
            p(paste("Showing", length(lvls), "color pickers (scroll to see all)")),
            color_pickers)
      } else {
        color_pickers
      }
    })
  })
  
  # --- Heatmap Logic ---
  observe({
    req(seurat_obj())
    obj <- seurat_obj()
    all_features <- c(rownames(obj), colnames(obj@meta.data))
    updateSelectizeInput(session, "hm_features", choices = all_features, server = TRUE)
    updateSelectInput(session, "hm_group_by", choices = names(obj@meta.data), selected = "seurat_clusters")
  })
  
  # Manual annotation color UI
  output$hm_manual_color_ui <- renderUI({
    req(seurat_obj(), input$hm_group_by)
    req(input$hm_anno_color_source == "Manual")
    
    obj <- seurat_obj()
    grp <- input$hm_group_by
    lvls <- sort(unique(as.character(obj@meta.data[[grp]])))
    
    if (length(lvls) > 50) return(p("Too many levels for manual colors"))
    
    lapply(lvls, function(l) {
      colourInput(paste0("hm_col_", l), paste("Color", l), value = "#808080")
    })
  })
  
  heatmap_obj <- reactiveVal(NULL)
  
  observeEvent(input$run_hm, {
    req(seurat_obj(), input$hm_features)
    
    withProgress(message = "Generating Heatmap...", value = 0, {
      tryCatch({
        obj <- seurat_obj()
        feats <- input$hm_features
        
        # Check which features exist in the object
        existing_feats <- feats[feats %in% c(rownames(obj), colnames(obj@meta.data))]
        
        # Get data matrix from 'data' layer
        # Handle features that might be metadata vs genes
        gene_feats <- existing_feats[existing_feats %in% rownames(obj)]
        meta_feats <- existing_feats[existing_feats %in% colnames(obj@meta.data)]
        
        if (length(gene_feats) == 0) {
          showNotification("Heatmap currently supports gene features from the 'data' layer", type="error")
          return()
        }
        
        # Extract data
        mat <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[gene_feats, , drop = FALSE])
        
        # Apply manual scaling
        scaling <- input$hm_scaling
        legend_title <- "Expression"
        if (scaling == "row") {
          # Row Z-score: (x - mean(row)) / sd(row)
          mat <- t(apply(mat, 1, function(x) {
            if(sd(x) == 0) return(x - mean(x))
            (x - mean(x)) / sd(x)
          }))
          legend_title <- "Row Z-score"
        }
        # Determine palette and mapping for viridis options
        pal <- input$hm_palette
        option_map <- list(viridis="D", magma="A", plasma="C", inferno="B", cividis="E")
        
        # Custom Heatmap using DoHeatmap requires a Seurat object with the data in scale.data
        # We create a minimal temporary object for this purpose
        tmp_obj <- subset(obj, features = gene_feats)
        tmp_obj <- SetAssayData(tmp_obj, layer = "scale.data", new.data = mat)
        
        # 2. Determine annotation colors FIRST
        anno_source <- input$hm_anno_color_source
        grp <- input$hm_group_by
        lvls <- sort(unique(as.character(obj@meta.data[[grp]])))
        anno_cols <- NULL
        
        if (anno_source == "Palette") {
          anno_pal <- input$hm_anno_palette
          n <- length(lvls)
          if (anno_pal %in% names(option_map)) {
            anno_cols <- viridis::viridis(n, option = option_map[[anno_pal]])
          } else if (anno_pal %in% names(MetBrewer::MetPalettes)) {
            anno_cols <- MetBrewer::met.brewer(anno_pal, n)
          } else {
            max_n <- RColorBrewer::brewer.pal.info[anno_pal, "maxcolors"]
            if (n <= max_n) {
              anno_cols <- RColorBrewer::brewer.pal(n, anno_pal)
            } else {
              anno_cols <- colorRampPalette(RColorBrewer::brewer.pal(max_n, anno_pal))(n)
            }
          }
          if (!is.null(anno_cols)) names(anno_cols) <- lvls
        } else if (anno_source == "Manual") {
          anno_cols <- sapply(lvls, function(l) {
            c <- input[[paste0("hm_col_", l)]]
            if(is.null(c)) "#808080" else c
          })
          names(anno_cols) <- lvls
        }

        # 3. Generate Heatmap with correct colors and gap size
        # DEFINITIVE FIX: Override standard segment defaults to white globally for this call
        orig_segment_color <- get_geom_defaults("segment")$colour
        update_geom_defaults("segment", list(colour = "white"))
        
        p <- tryCatch({
          DoHeatmap(
            tmp_obj, 
            features = gene_feats, 
            group.by = input$hm_group_by, 
            group.colors = anno_cols, 
            draw.lines = TRUE, 
            lines.width = input$hm_lines_width, 
            raster = input$hm_raster,
            size = input$hm_font_size / 2.5
          )
        }, finally = {
          # Restore original default
          update_geom_defaults("segment", list(colour = orig_segment_color))
        })
        
        if (!is.null(p)) {
          # 4. Apply expression color palette
          if (pal %in% names(option_map)) {
            p <- p + scale_fill_viridis(option = option_map[[pal]], name = legend_title)
          } else {
            p <- p + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = pal)), name = legend_title)
          }
          
          # 5. Fix: Change separation lines from grey to white and set background to white
          # Fallback: Intensive search through all list elements for grey/gray and converting to white
          p$layers <- lapply(p$layers, function(l) {
            # 1. Check aes_params
            if (!is.null(l$aes_params$colour) && is.character(l$aes_params$colour)) {
              if (tolower(l$aes_params$colour) %in% c("grey", "gray", "grey50", "gray50", "#7f7f7f", "#808080", "grey80", "gray80")) {
                l$aes_params$colour <- "white"
              }
            }
            # 2. Check geom_params if it exists
            if (!is.null(l$geom_params$colour) && is.character(l$geom_params$colour)) {
              if (tolower(l$geom_params$colour) %in% c("grey", "gray", "grey50", "gray50", "#7f7f7f", "#808080")) {
                l$geom_params$colour <- "white"
              }
            }
            return(l)
          })
          
          p <- p + theme(
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.grid = element_blank(),
            # Fonts
            axis.text.y = element_text(size = input$hm_font_size),
            legend.text = element_text(size = input$hm_legend_font_size),
            legend.title = element_text(size = input$hm_legend_font_size)
          )
        }
        
        heatmap_obj(p)
        showNotification("Heatmap generated!", type="message")
      }, error = function(e) {
        showNotification(paste("Heatmap error:", e$message), type="error")
      })
    })
  })
  
  output$heatmap_plot <- renderPlot({
    req(heatmap_obj())
    heatmap_obj()
  })
  
  output$download_hm <- downloadHandler(
    filename = function() { paste0("heatmap_", Sys.Date(), ".pdf") },
    content = function(file) {
      req(heatmap_obj())
      ggsave(file, heatmap_obj(), width = input$hm_width, height = input$hm_height, units = "in")
    }
  )
  
  # --- Render Plots ---
  plot_list <- reactiveValues()
  
  observe({
    req(seurat_obj())
    plot_list$p1 <- tryCatch(generate_plot("p1", seurat_obj()), 
                             error=function(e) ggplot() + annotate("text", x=0, y=0, label="Configure plot settings", size=5, color="gray50") + theme_void())
  })
  observe({
    req(seurat_obj())
    plot_list$p2 <- tryCatch(generate_plot("p2", seurat_obj()), 
                             error=function(e) ggplot() + annotate("text", x=0, y=0, label="Configure plot settings", size=5, color="gray50") + theme_void())
  })
  observe({
    req(seurat_obj())
    plot_list$p3 <- tryCatch(generate_plot("p3", seurat_obj()), 
                             error=function(e) ggplot() + annotate("text", x=0, y=0, label="Configure plot settings", size=5, color="gray50") + theme_void())
  })
  observe({
    req(seurat_obj())
    plot_list$p4 <- tryCatch(generate_plot("p4", seurat_obj()), 
                             error=function(e) ggplot() + annotate("text", x=0, y=0, label="Configure plot settings", size=5, color="gray50") + theme_void())
  })
  
  output$plot1 <- renderPlot(plot_list$p1)
  output$plot2 <- renderPlot(plot_list$p2)
  output$plot3 <- renderPlot(plot_list$p3)
  output$plot4 <- renderPlot(plot_list$p4)
  
  # --- Export ---
  output$download_plot <- downloadHandler(
    filename = function() paste0("seurat_plot_", Sys.Date(), ".", input$export_format),
    content = function(file) {
      tgt <- input$export_target
      p <- if(tgt=="All Plots") (plot_list$p1 + plot_list$p2) / (plot_list$p3 + plot_list$p4)
           else if(tgt=="Plot 1") plot_list$p1
           else if(tgt=="Plot 2") plot_list$p2
           else if(tgt=="Plot 3") plot_list$p3
           else if(tgt=="Plot 4") plot_list$p4
      
      ggsave(file, p, width=input$export_w, height=input$export_h, device=input$export_format)
    }
  )
}

shinyApp(ui, server)
