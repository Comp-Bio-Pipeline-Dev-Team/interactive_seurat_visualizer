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
library(shinyjs)

# Load NMF BEFORE BiocGenerics packages to avoid do.call masking
library(NMF)

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
    
    /* Dark Mode Styles */
    body.dark-mode { background-color: #2c3e50; color: #ecf0f1; }
    body.dark-mode .plot-container { background-color: #34495e; border-color: #7f8c8d; }
    body.dark-mode .well { background-color: #34495e; border: none; color: #ecf0f1; }
    body.dark-mode .nav-tabs > li > a { color: #ecf0f1; }
    body.dark-mode .nav-tabs > li.active > a, 
    body.dark-mode .nav-tabs > li.active > a:focus, 
    body.dark-mode .nav-tabs > li.active > a:hover { color: #2c3e50; background-color: #ecf0f1; }
      /* Dark Mode Table Fixes */
      body.dark-mode .dataTables_wrapper .dataTables_length, 
      body.dark-mode .dataTables_wrapper .dataTables_filter, 
      body.dark-mode .dataTables_wrapper .dataTables_info, 
      body.dark-mode .dataTables_wrapper .dataTables_processing, 
      body.dark-mode .dataTables_wrapper .dataTables_paginate {
        color: #e0e0e0 !important;
      }
      body.dark-mode table.dataTable tbody tr {
        background-color: #2b2b2b !important;
        color: #e0e0e0 !important;
      }
      body.dark-mode .dataTables_wrapper .dataTables_paginate .paginate_button {
        color: #e0e0e0 !important; 
      }
      
      /* Tabs */
      .nav-pills > li.active > a, .nav-pills > li.active > a:focus, .nav-pills > li.active > a:hover {
        color: #fff;
        background-color: #337ab7;
      }
    body.dark-mode h4, body.dark-mode h5 { color: #ecf0f1; }
  ")))
  ,tags$script(HTML('
    $(document).on("click", ".plot-container", function() {
      // 1. Highlight the container
      $(".plot-container").removeClass("active-plot");
      $(this).addClass("active-plot");
      
      // 2. Identify which plot it is
      var id = $(this).attr("id");
      var targetTab = "";
      if(id === "plot1-container") targetTab = "Plot 1";
      if(id === "plot2-container") targetTab = "Plot 2";
      if(id === "plot3-container") targetTab = "Plot 3";
      if(id === "plot4-container") targetTab = "Plot 4";
      
      // 3. Switch the settings tab in the sidebar
      // We find the tab link with the matching data-value and trigger a click
      if(targetTab !== "") {
        var selector = "a[data-value=\'" + targetTab + "\']";
        $(selector).tab("show");
      }
    });
  ')),
  
  # Dark Mode Toggle in Footer or Header
  footer = div(
    style = "position: fixed; bottom: 10px; right: 10px; z-index: 1000;",
    checkboxInput("dark_mode", "Dark Mode", value = FALSE)
  ),
  
  
  
  # Initialize shinyjs
  useShinyjs(),
  
  tabPanel("Quality Control",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("1. Load Data"),
        fileInput("seurat_file", "Upload Seurat (.rds or .h5ad)", accept = c(".rds", ".h5ad")),
        radioButtons("ensembl_species", "Species", choices = c("Human", "Mouse"), inline = TRUE),
        actionButton("convert_ensembl", "Convert Ensembl IDs to Symbols", icon = icon("dna"), class = "btn-info btn-block"),
        hr(),
        h4("2. Filtering"),
        h5("QC Thresholds"),
        selectInput("qc_mito_col", "Mitochondrial Column", choices = NULL),
        sliderInput("filter_percent_mt", "Max % Mito", min=0, max=100, value=20),
        sliderInput("filter_nFeature_min", "Min nFeature", min=0, max=10000, value=200),
        sliderInput("filter_nCount_min", "Min nCount", min=0, max=50000, value=500),
        hr(),
        h5("Metadata Subset"),
        selectInput("subset_col", "Meta Col", choices = NULL),
        selectizeInput("subset_levels", "Keep Levels", choices = NULL, multiple=TRUE),
        br(),
        actionButton("apply_filters", "Apply Filtering", class="btn-warning btn-block"),
        br(),
        actionButton("reset_data", "Reset to Original", class="btn-danger btn-block"),
        hr(),
        h5("Visualization"),
        p("Distribution of QC metrics across samples/groups.")
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel("Violin Plots", plotOutput("qc_violin", height = "600px")),
          tabPanel("Scatter Plot", plotOutput("qc_scatter", height = "600px"))
        )
      )
    )
  ),
  
  tabPanel("Signatures",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("UCell Signatures"),
        p("Calculate signature scores for gene sets."),
        textAreaInput("sig_genes", "Genes (one per line)", rows=10, placeholder="CD3D\nCD3E\nCD8A"),
        textInput("sig_name", "Signature Name", "MySig"),
        actionButton("calc_ucell", "Calculate UCell", class="btn-success btn-block")
      ),
      mainPanel(
         width = 9,
         h4("Instructions"),
         p("1. Enter a list of genes."),
         p("2. Give the signature a name."),
         p("3. Click Calculate."),
         p("4. The score will be added as metadata. Go to 'Visualization' to plot it."),
         hr(),
         verbatimTextOutput("ucell_log")
      )
    )
  ),
  
  tabPanel("Visualization",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        p(tags$small("Data loaded in 'Quality Control' tab.")),
        br(),
        tabsetPanel(id = "settings_tabs",
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
          column(6, div(style="text-align: center; font-weight: bold; color: #2c3e50;", "Plot 1"), div(id="plot1-container", class="plot-container", plotOutput("plot1", height="400px"))),
          column(6, div(style="text-align: center; font-weight: bold; color: #2c3e50;", "Plot 2"), div(id="plot2-container", class="plot-container", plotOutput("plot2", height="400px")))
        ),
        fluidRow(
          column(6, div(style="text-align: center; font-weight: bold; color: #2c3e50;", "Plot 3"), div(id="plot3-container", class="plot-container", plotOutput("plot3", height="400px"))),
          column(6, div(style="text-align: center; font-weight: bold; color: #2c3e50;", "Plot 4"), div(id="plot4-container", class="plot-container", plotOutput("plot4", height="400px")))
        )
      )
    )
  ),

  tabPanel("Gene Programs",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("NMF Settings"),
        sliderInput("nmf_k", "Number of Factors (k)", min=2, max=20, value=5),
        numericInput("nmf_nrun", "N Runs", value=10, min=1),
        actionButton("run_nmf", "Run NMF", class="btn-primary btn-block"),
        hr(),
        conditionalPanel(
          condition = "output.nmf_done",
          h4("Visualization"),
          selectInput("nmf_viz_type", "Plot Type", choices=c("Factor Heatmap", "FeaturePlot (Scores)")),
          selectInput("nmf_factor", "Select Factor", choices=NULL),
          sliderInput("nmf_pt_size", "Point Size (FeaturePlot)", min=0.1, max=3, value=1, step=0.1),
          selectInput("nmf_color", "Color Scale", choices=c("Magma", "Viridis", "Plasma", "Inferno"), selected="Magma"),
          selectInput("nmf_reduction", "Reduction", choices = c("umap", "tsne", "pca"), selected = "umap"),
          downloadButton("download_nmf_csv", "Download Top Genes")
        )
      ),
      mainPanel(
        width = 9,
        conditionalPanel(
          condition = "input.nmf_viz_type == 'Factor Heatmap'",
          plotOutput("nmf_heatmap", height="600px")
        ),
        conditionalPanel(
          condition = "input.nmf_viz_type == 'FeaturePlot (Scores)'",
          plotOutput("nmf_featureplot", height="600px")
        ),
        hr(),
        h5("Top Genes for Selected Factor"),
        DT::dataTableOutput("nmf_gene_table")
      )
    )
  ),

  tabPanel("Differential Expression",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Differential Expression"),
        selectInput("de_group_by", "Comparison Group", choices = NULL),
        selectInput("de_ident_1", "Ident 1", choices = NULL),
        selectInput("de_ident_2", "Ident 2", choices = NULL),
        selectInput("de_test", "Test Method", choices = c("wilcox", "bimod", "t", "roc"), selected = "wilcox"),
        numericInput("de_logfc", "Min LogFC", value = 0.25, step = 0.05),
        numericInput("de_minpct", "Min Pct", value = 0.1, step = 0.05),
        numericInput("de_pval", "P-value Threshold", value = 0.05, step = 0.01),
        actionButton("run_de", "Run Differential Expression", class = "btn-primary btn-block"),
        hr(),
        h5("Volcano Plot Settings"),
        numericInput("volcano_pt_size", "Point Size", value = 1.5, min = 0.1),
        numericInput("volcano_alpha", "Alpha", value = 0.6, min = 0.1, max = 1),
        colourpicker::colourInput("volcano_sig_color", "Significant Color", value = "#E74C3C"),
        colourpicker::colourInput("volcano_nonsig_color", "Non-significant Color", value = "#95A5A6"),
        selectInput("volcano_format", "Export Format", choices = c("png", "pdf"), selected = "png"),
        numericInput("volcano_width", "Width", 8),
        numericInput("volcano_height", "Height", 6),
        downloadButton("download_volcano", "Download Volcano", class = "btn-info btn-block"),
        hr(),
        downloadButton("download_de", "Download Results Table", class = "btn-info btn-block")
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel("Results Table", DT::DTOutput("de_table")),
          tabPanel("Volcano Plot", plotlyOutput("de_volcano", height = "600px"))
        )
      )
    )
  ),

  tabPanel("Pathway Enrichment",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Enrichment Analysis"),
        selectInput("enrich_organism", "Organism", choices = c("Human" = "human", "Mouse" = "mouse"), selected = "human"),
        selectInput("enrich_source", "Gene Source", 
                    choices = c("From DE Results" = "de", "Custom List" = "custom")),
        
        conditionalPanel(
          condition = "input.enrich_source == 'de'",
          selectInput("enrich_gene_direction", "Direction", 
                      choices = c("All Significant" = "all", "Upregulated" = "up", "Downregulated" = "down"))
        ),
        
        conditionalPanel(
          condition = "input.enrich_source == 'custom'",
          textAreaInput("enrich_custom_genes", "Paste Genes (Symbol)", rows = 5, placeholder = "TP53\nEGFR\nCD4")
        ),
        
        selectInput("enrich_method", "Method", 
                    choices = c("Over-Representation (ORA)" = "ora", "GSEA" = "gsea")),
        
        selectInput("enrich_db", "Database", 
                    choices = c("GO Biological Process" = "go_bp", 
                                "GO Molecular Function" = "go_mf",
                                "KEGG Pathways" = "kegg", 
                                "Reactome" = "reactome",
                                "MSigDB Hallmark" = "hallmark")),
        
        actionButton("run_enrichment", "Run Enrichment Analysis", class = "btn-success btn-block"),
        hr(),
        h5("Visualization"),
        selectInput("enrich_plot_type", "Plot Type", 
                    choices = c("Dotplot", "Barplot", "Network (Cnet)", "Enrichment Map (Emap)", "GSEA Plot")),
        numericInput("enrich_show_n", "Show Top N", value = 15, min = 1, max = 50),
        downloadButton("download_enrich_plot", "Download Plot")
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel("Enrichment Plot", plotOutput("enrich_plot", height = "700px")),
          tabPanel("Results Table", DT::DTOutput("enrich_table"))
        )
      )
    )
  ),

  tabPanel("Heatmap",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Heatmap Settings"),
        selectizeInput("hm_features", "Features (Genes)", choices = NULL, multiple = TRUE),
        selectInput("hm_group_by", "Group By", choices = NULL),
        checkboxInput("hm_cluster_rows", "Cluster Rows", value = TRUE),
        checkboxInput("hm_cluster_cols", "Cluster Columns", value = FALSE),
        selectInput("hm_scaling", "Scaling", choices = c("Row Z-Score" = "row", "None" = "none"), selected = "row"),
        selectInput("hm_palette", "Palette", choices = c("Viridis", "Magma", "RdBu", "RdYlBu"), selected = "Viridis"),
        
        hr(),
        h5("Annotation Colors"),
        selectInput("hm_anno_color_source", "Color Source", choices = c("Default", "Manual"), selected = "Default"),
        conditionalPanel(
          condition = "input.hm_anno_color_source == 'Manual'",
          uiOutput("hm_manual_color_ui")
        ),
        
        hr(),
        actionButton("generate_heatmap", "Generate Heatmap", class = "btn-primary btn-block"),
        br(),
        downloadButton("download_heatmap", "Download Heatmap", class = "btn-info btn-block")
      ),
      mainPanel(
        width = 9,
        plotOutput("heatmap_plot", height = "800px")
      )
    )
  ),

tabPanel("Multimodal (ADT)",
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Protein Expression"),
      selectizeInput("adt_feature", "Select Protein", choices = NULL),
      sliderInput("pt_size_adt", "Point Size", min=0.1, max=3, value=1),
      selectInput("adt_palette", "Color Scale", choices = c("Viridis" = "viridis", "Magma" = "magma", "Plasma" = "plasma", "Inferno" = "inferno", "Cividis" = "cividis", "Spectral" = "Spectral", "RdBu" = "RdBu", "RdYlBu" = "RdYlBu"), selected = "viridis"),
      selectInput("species_select", "Species", choices = c("Human", "Mouse"), selected = "Human"),
      hr(),
      h4("Co-expression"),
      selectizeInput("adt_feature_x", "Protein (X-axis)", choices = NULL),
      selectizeInput("rna_feature_y", "Gene (Y-axis)", choices = NULL),
      hr(),
      h5("Dimensions & Aesthetics"),
      fluidRow(
        column(6, numericInput("adt_export_width", "Width (in)", value = 8, min = 1, max = 50)),
        column(6, numericInput("adt_export_height", "Height (in)", value = 6, min = 1, max = 50))
      ),
      selectInput("adt_export_format", "Format", choices = c("png", "pdf", "jpg"), selected = "png"),
      numericInput("adt_pt_size", "Point Size", min=0.1, max=5, value=1, step = 0.1),
      numericInput("adt_font_size", "Font Size", value = 12, min = 1, max = 24),
      br(),
      downloadButton("download_adt_feature", "Download FeaturePlot", class="btn-info btn-block"),
      downloadButton("download_adt_coexp", "Download Co-expression", class="btn-info btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("FeaturePlot", plotOutput("adt_feature_plot", height = "600px")),
        tabPanel("Co-expression", plotOutput("adt_rna_coexpression", height = "600px"))
      )
    )
  )
),

tabPanel("Gating",
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Co-expression Gating"),
      p("Select two features (e.g., ADT proteins) and gate cells."),
      selectizeInput("gate_feature_x", "Feature X (e.g. CD3)", choices = NULL),
      selectizeInput("gate_feature_y", "Feature Y (e.g. CD19)", choices = NULL),
      hr(),
      h5("Gating Actions"),
      p(tags$small("1. Use Lasso/Box tool on plot to select cells.")),
      p(tags$small("2. Click 'Gate' to filter dataset.")),
      actionButton("gate_cells", "Gate Selected Cells", class="btn-danger btn-block", icon=icon("filter")),
      br(),
      actionButton("reset_gating", "Reset Gating", class="btn-secondary btn-block")
    ),
    mainPanel(
      width = 9,
      plotlyOutput("gating_plot", height = "600px"),
      verbatimTextOutput("gating_info")
    )
  )
),

tabPanel("VDJ Analysis",
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Upload VDJ Data"),
      fileInput("vdj_upload", "Upload contigs.csv", accept = c(".csv")),
      p(tags$small("Upload 'filtered_contig_annotations.csv' from CellRanger.")),
      hr(),
      h4("Visualization"),
      selectInput("vdj_viz_type", "Plot Type", 
                  choices = c("Clonal Overlay", "V Gene Usage", "Clonal Diversity"), 
                  selected = "Clonal Overlay"),
      
      conditionalPanel(
        condition = "input.vdj_viz_type == 'Clonal Diversity'",
        selectInput("vdj_group_by", "Group By", choices = NULL),
        selectInput("vdj_diversity_metric", "Metric", 
                    choices = c("Shannon", "Simpson", "inv.simpson", "chao"), 
                    selected = "Shannon"),
        numericInput("vdj_boot_n", "Bootstraps", value = 100, min = 10, max = 1000)
      ),
      hr(),
      h5("Export"),
      downloadButton("download_vdj_plot", "Download Plot", class="btn-info btn-block")
    ),
    mainPanel(
      width = 9,
      conditionalPanel(
        condition = "input.vdj_viz_type == 'Clonal Overlay'",
        plotOutput("vdj_overlay", height = "600px")
      ),
      conditionalPanel(
        condition = "input.vdj_viz_type == 'V Gene Usage'",
        plotOutput("vdj_vgene", height = "600px")
      ),
      conditionalPanel(
        condition = "input.vdj_viz_type == 'Clonal Diversity'",
        plotOutput("vdj_diversity_plot", height = "600px")
      )
    )
  )
),

tabPanel("Pseudobulk",
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Pseudobulk Settings"),
      selectInput("pb_group_by", "Group By", choices = NULL),
      p(tags$small("Aggregates expression by the selected metadata column.")),
      checkboxInput("pb_cluster_rows", "Cluster Rows", value = TRUE),
      checkboxInput("pb_cluster_cols", "Cluster Cols", value = TRUE),
      selectInput("pb_scaling", "Scaling", choices = c("Row Z-Score"="row", "None"="none"), selected="row"),
      selectInput("pb_palette", "Palette", choices = c("Viridis", "Magma", "RdBu", "RdYlBu"), selected="Viridis"),
      hr(),
      h5("Export"),
      fluidRow(
          column(6, numericInput("pb_export_w", "W", 10)),
          column(6, numericInput("pb_export_h", "H", 10))
      ),
      downloadButton("download_pb_heatmap", "Download Heatmap", class="btn-info btn-block"),
      hr(),
      p("This visualization helps identify robust group-level differences.")
    ),
    mainPanel(
      width = 9,
      plotOutput("pseudobulk_heatmap", height = "800px")
    )
  )
)
)

server <- function(input, output, session) {
  
  # Initialize Reactives
  original_obj <- reactiveVal(NULL)
  # Main reactive holding the (potentially subsetted) Seurat object
  seurat_obj <- reactiveVal(NULL)
  
  # Alias for backend scripts which expect 'data_object()'
  data_object <- seurat_obj

  # Source Backend Logic (Local Scope)
  source("multimodal_backend.R", local = TRUE)
  source("vdj_backend.R", local = TRUE)
  source("qc_backend.R", local = TRUE)
  source("pseudobulk_backend.R", local = TRUE)
  source("enrichment_backend.R", local = TRUE)
  
  # Dark Mode Observer
  observe({
    if (input$dark_mode) {
      runjs("$('body').addClass('dark-mode');")
    } else {
      runjs("$('body').removeClass('dark-mode');")
    }
  })
  
  observeEvent(input$seurat_file, {
    req(input$seurat_file)
    tryCatch({
      path <- input$seurat_file$datapath
      ext <- tolower(tools::file_ext(input$seurat_file$name))
      
      obj <- NULL
      if (ext == "rds") {
        obj <- readRDS(path)
      } else if (ext == "h5ad") {
        # Strategy: Zellkonverter (User Preferred)
        if (requireNamespace("zellkonverter", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
          showNotification("Reading h5ad with zellkonverter...", type="message", duration=2)
          tryCatch({
            sce <- zellkonverter::readH5AD(path)
            obj <- as.Seurat(sce, counts = "X", data = "X")
          }, error = function(e) {
            stop(paste("Zellkonverter loading failed:", e$message))
          })
        } else {
          stop("zellkonverter package is missing. Please install it to load .h5ad files.")
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
      
      # Safety Check: Run NormalizeData if max value suggests raw counts > 20
      # Use tryCatch to avoid crashing if memory is tight
      tryCatch({
        data_mat <- GetAssayData(obj, layer = "data")
        if (max(data_mat, na.rm=TRUE) > 20) {
          showNotification("Detected raw counts. Normalizing data...", type = "message", duration = 5)
          obj <- NormalizeData(obj)
        }
      }, error = function(e) {
        warning("Normalization check failed: ", e$message)
      })
      
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
      
      # Calculate QC metrics if missing
      if (!"percent.mt" %in% colnames(obj@meta.data)) {
        if (any(grepl("^MT-", rownames(obj)))) {
            obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
        } else if (any(grepl("^mt-", rownames(obj)))) {
            obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
        } else {
            obj[["percent.mt"]] <- 0
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
    updateSelectInput(session, "pb_group_by", choices = meta_cols, selected = "seurat_clusters") # For Pseudobulk
    
    # Update Mitochondrial Column choices (numeric only)
    numeric_cols <- names(seurat_obj()@meta.data)[sapply(seurat_obj()@meta.data, is.numeric)]
    
    # Check for likely candidates
    selected_mito <- "percent.mt"
    if (!selected_mito %in% numeric_cols) {
       candidates <- grep("^percent\\.|^pct_|^mt", numeric_cols, value = TRUE, ignore.case = TRUE)
       if(length(candidates) > 0) selected_mito <- candidates[1]
       else selected_mito <- numeric_cols[1] # Default to first numeric
    }
    updateSelectInput(session, "qc_mito_col", choices = numeric_cols, selected = selected_mito)
  })

  # Load Demo Data Observer
  observeEvent(input$load_demo, {
    req(file.exists("test_seurat.rds"))
    
    withProgress(message = 'Loading Demo Data...', value = 0, {
      tryCatch({
        obj <- readRDS("test_seurat.rds")
        
        # Ensure percent.mt
        if (!"percent.mt" %in% colnames(obj@meta.data)) {
           if (any(grepl("^MT-", rownames(obj)))) {
             obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
           } else if (any(grepl("^mt-", rownames(obj)))) {
             obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
           } else {
             obj[["percent.mt"]] <- 0
           }
        }
        
        # Update reactive values
        seurat_obj(obj)
        original_obj(obj)
        
        # Update UI choices
        meta_cols <- names(obj@meta.data)
        updateSelectInput(session, "subset_col", choices = meta_cols)
        
        # Populate Gating features (Prefer ADT if available)
        features <- rownames(obj)
        if ("ADT" %in% names(obj@assays)) {
             adt_feats <- rownames(obj[["ADT"]])
             # Prefix ADT features for clarity if needed, or just list them
             features <- c(adt_feats, features)
        }
        updateSelectizeInput(session, "gate_feature_x", choices = features, server = TRUE)
        updateSelectizeInput(session, "gate_feature_y", choices = features, server = TRUE)
        
        showNotification("Demo data loaded successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Error loading demo:", e$message), type = "error")
      })
    })
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
           subtitle = paste("Positive LFC = Upregulated in", input$de_ident_1),
           x = "avg_log2FC", y = "-log10(p-value)")
  })
  
  output$de_volcano <- renderPlotly({
    req(volcano_plot())
    ggplotly(volcano_plot(), tooltip = "text")
  })
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      # Sanitize filename components to prevent download failures
      i1 <- gsub("[^[:alnum:]]", "_", input$de_ident_1)
      i2 <- if(is.null(input$de_ident_2)) "Others" else gsub("[^[:alnum:]]", "_", input$de_ident_2)
      paste0("Volcano_", i1, "_vs_", i2, "_", Sys.Date(), ".", input$volcano_format) 
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
  
  observeEvent(input$apply_filters, {
    req(original_obj())
    
    # Start with original object to allow re-filtering
    obj <- original_obj()
    
    # 1. Apply Numeric QC Filters
    # Check if percent.mt exists
    if (!"percent.mt" %in% colnames(obj@meta.data)) {
      if (any(grepl("^MT-", rownames(obj)))) {
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
      } else {
        # If no MT genes, create dummy 0 column to avoid error if slider is used (or warn)
        obj[["percent.mt"]] <- 0
        showNotification("No mitochondrial genes found (pattern ^MT-). percent.mt set to 0.", type="warning")
      }
    }
    
    # Check nFeature/nCount existence (standard Seurat)
    # Usually nFeature_RNA, nCount_RNA. If integrated, might vary.
    # We'll use the default assay's meta features if possible
    meta <- obj@meta.data
    
    # Filter logic
    keep_cells <- colnames(obj)
    
    # Filter by percent.mt (User selected column)
    mito_col <- input$qc_mito_col
    if (!is.null(mito_col) && mito_col %in% colnames(meta)) {
       cells_mt <- rownames(meta)[meta[[mito_col]] <= input$filter_percent_mt]
       keep_cells <- intersect(keep_cells, cells_mt)
    }
    
    # Filter by nFeature
    nf_col <- paste0("nFeature_", DefaultAssay(obj))
    if (!nf_col %in% colnames(meta)) nf_col <- "nFeature_RNA" # Fallback
    if (nf_col %in% colnames(meta)) {
       cells_nf <- rownames(meta)[meta[[nf_col]] >= input$filter_nFeature_min]
       keep_cells <- intersect(keep_cells, cells_nf)
    }
    
    # Filter by nCount
    nc_col <- paste0("nCount_", DefaultAssay(obj))
    if (!nc_col %in% colnames(meta)) nc_col <- "nCount_RNA" # Fallback
    if (nc_col %in% colnames(meta)) {
       cells_nc <- rownames(meta)[meta[[nc_col]] >= input$filter_nCount_min]
       keep_cells <- intersect(keep_cells, cells_nc)
    }
    
    # 2. Apply Metadata Subset (if selected)
    if (!is.null(input$subset_col) && !is.null(input$subset_levels)) {
       cells_meta <- rownames(meta)[meta[[input$subset_col]] %in% input$subset_levels]
       keep_cells <- intersect(keep_cells, cells_meta)
    }
    
    # Apply subset
    if (length(keep_cells) == 0) {
      showNotification("Filtering resulted in 0 cells! Reverting...", type="error")
    } else {
      # Use subset function
      obj_sub <- subset(obj, cells = keep_cells)
      seurat_obj(obj_sub)
      showNotification(paste("Filtered to", length(keep_cells), "cells"), type="message")
    }
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
  
  # --- Color Helper ---
  get_colors <- function(id, obj, group_var, for_scpubr = FALSE) {
    ns <- function(x) paste0(id, "-", x)
    source <- input[[ns("color_source")]]
    
    if (source == "Default") return(NULL)
    
    # Determine levels
    if (is.null(group_var) || group_var == "None" || group_var == "Default") {
      lvls <- levels(Idents(obj))
    } else {
      lvls <- sort(unique(obj@meta.data[[group_var]]))
    }
    
    n <- length(lvls)
    if (n == 0) return(NULL)
    
    if (source == "Palette") {
      pal <- input[[ns("palette_name")]]
      
      # Viridis family - fix option mapping
      if (pal %in% c("viridis", "magma", "plasma", "inferno", "cividis")) {
        option_map <- list(viridis="D", magma="A", plasma="C", inferno="B", cividis="E")
        cols <- viridis::viridis(n, option=option_map[[pal]])
      } else if (pal %in% names(MetBrewer::MetPalettes)) {
        # MetBrewer
        cols <- MetBrewer::met.brewer(pal, n)
      } else if (pal %in% c("Set1", "Set2", "Set3", "Dark2", "Paired", "Pastel1", "Accent")) {
        # RColorBrewer
        max_n <- RColorBrewer::brewer.pal.info[pal, "maxcolors"]
        if (n <= max_n) {
          cols <- RColorBrewer::brewer.pal(n, pal)
        } else {
          cols <- colorRampPalette(RColorBrewer::brewer.pal(max_n, pal))(n)
        }
      } else {
        return(NULL)
      }
      
      # For SCpubr, return named vector
      if (for_scpubr) {
        names(cols) <- lvls
      }
      return(cols)
    }
    
    if (source == "Manual") {
      cols <- sapply(lvls, function(l) {
        c <- input[[ns(paste0("col_", l))]]
        if(is.null(c)) "gray" else c
      })
      # For SCpubr, ensure it's a named vector
      if (for_scpubr) {
        names(cols) <- lvls
      }
      return(cols)
    }
    
    return(NULL)
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
    
    colors <- get_colors(id, obj, if(is.null(grp)) "Default" else grp)
    # For SCpubr, get named colors
    style <- input[[ns("plot_style")]]
    if (!is.null(style) && style == "SCpubr") {
      colors <- get_colors(id, obj, if(is.null(grp)) "Default" else grp, for_scpubr = TRUE)
    }
    pt_size <- input[[ns("pt_size")]]
    if (is.null(pt_size)) pt_size <- 1  # Default if not set
    
    p <- NULL
    
    if (ptype == "ClusterDistrBar") {
      g1 <- input[[ns("cdb_group1")]]
      g2 <- input[[ns("cdb_group2")]]
      req(g1, g2)
      
      # Get colors for the fill variable (g2), not the x-axis
      colors <- get_colors(id, obj, g2)
      
      df <- as.data.frame(table(obj@meta.data[[g1]], obj@meta.data[[g2]]))
      colnames(df) <- c("Sample", "Cluster", "Count")
      
      p <- ggplot(df, aes(x=Sample, y=Count, fill=Cluster)) + 
        geom_bar(stat="identity", position="fill") +
        scale_y_continuous(labels = scales::percent) +
        labs(y="Proportion", x=g1, fill=g2) +
        theme_classic() +
        theme(axis.text.x = element_text(angle=45, hjust=1))
      
      if (!is.null(colors)) p <- p + scale_fill_manual(values=colors)
      
      # Add cell counts if requested
      show_counts <- input[[ns("show_counts")]]
      if (!is.null(show_counts) && show_counts) {
        # Calculate totals per sample
        df_totals <- aggregate(Count ~ Sample, df, sum)
        colnames(df_totals) <- c("Sample", "Total")
        df_totals$label_y <- 1.05  # Position above bar
        
        # Get count size (default to 3 if not set)
        count_size <- input[[ns("count_size")]]
        if (is.null(count_size)) count_size <- 3
        
        p <- p + geom_text(data=df_totals, aes(x=Sample, y=label_y, label=paste0("n=", Total)), 
                          inherit.aes=FALSE, size=count_size, fontface="bold")
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
            # FeaturePlot - Enforce viridis for sequential data as requested
            # We ignore other palettes/manual colors for sequential data in SCpubr mode
            
            # Default to 'viridis' (D)
            viridis_pal <- "D"
            viridis_dir <- 1 
            
            # Check if user selected a specific viridis option
            pal_source <- input[[ns("color_source")]]
            if (!is.null(pal_source) && pal_source == "Palette") {
              pal_name <- input[[ns("palette_name")]]
              # Respect Viridis options
              if (!is.null(pal_name) && pal_name %in% c("viridis", "magma", "plasma", "inferno", "cividis")) {
                 option_map <- list(viridis="D", magma="A", plasma="C", inferno="B", cividis="E")
                 viridis_pal <- option_map[[pal_name]]
              }
              # If user selected RColorBrewer/MetBrewer, we can't easily force it into SCpubr's 'viridis_color_map' argument
              # which typically expects a viridis option character. 
              # However, SCpubr might support 'colors.use' for continuous? 
              # Actually SCpubr do_FeaturePlot uses 'viridis_color_map' = character.
              # If we strictly want to support others, we might need to fallback to standard, 
              # OR just warn user. For now, we stick to respecting viridis variants.
            }

            
            # Get direction
            if (!is.null(input[[ns("scpubr_viridis_dir")]])) {
               viridis_dir <- as.numeric(input[[ns("scpubr_viridis_dir")]])
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
              use_viridis = TRUE,
              viridis.palette = viridis_pal,
              viridis.direction = viridis_dir
            )
          } else if (ptype == "ViolinPlot") {
            req(feat)
             # ViolinPlot is tricky: can be discrete (if grouped) or just density.
             # SCpubr::do_ViolinPlot usually uses discrete colors for groups.
             # We pass 'colors.use' which handles all discrete palettes (Viridis/Brewer/Manual)
             
            p <- tryCatch({
              SCpubr::do_ViolinPlot(
                sample = obj,
                features = feat,
                group.by = grp,
                split.by = splt,
                colors.use = colors, # This supports the full range of discrete palettes
                font.size = font_size,
                legend.position = legend_pos
              )
            }, error = function(e) {
               # If it still fails (e.g. package issues), show error but cleaner
               ggplot() + 
                 annotate("text", x=0.5, y=0.5, label=paste("SCpubr Violin Error:\n", e$message), size=5, color="red") + 
                 theme_void() + 
                 theme(plot.margin = margin(10,10,10,10))
            })
          } else if (ptype == "DotPlot") {
            req(feat)
            # DotPlot is sequential (expression) + discrete (size). 
            # Enforce viridis for the color scale.
            
            viridis_pal <- "E" # Default to cividis or viridis
            viridis_dir <- 1
            
            pal_source <- input[[ns("color_source")]]
            if (!is.null(pal_source) && pal_source == "Palette") {
              pal_name <- input[[ns("palette_name")]]
              if (!is.null(pal_name) && pal_name %in% c("viridis", "magma", "plasma", "inferno", "cividis")) {
                 viridis_map <- list(magma="A", inferno="B", plasma="C", viridis="D", cividis="E")
                 viridis_pal <- viridis_map[[pal_name]]
              }
            }
            
            if (!is.null(input[[ns("scpubr_viridis_dir")]])) {
               viridis_dir <- as.numeric(input[[ns("scpubr_viridis_dir")]])
            }
            
            p <- tryCatch({
              SCpubr::do_DotPlot(
                sample = obj,
                features = feat,
                group.by = grp,
                split.by = splt,
                font.size = font_size,
                legend.position = legend_pos,
                use_viridis = TRUE,
                viridis.palette = viridis_pal,
                viridis.direction = viridis_dir
              )
            }, error = function(e) {
              # If split.by fails, try without it
              result <- SCpubr::do_DotPlot(
                sample = obj,
                features = feat,
                group.by = grp,
                font.size = font_size,
                legend.position = legend_pos,
                use_viridis = TRUE,
                viridis.palette = viridis_pal,
                viridis.direction = viridis_dir
              )
            })
          }
          
          # Fix SCpubr Legend Orientation
          if (!is.null(p)) {
            leg_dir <- if (legend_pos %in% c("bottom", "top")) "horizontal" else "vertical"
            
            # Define guides based on plot type
            my_guides <- NULL
            if (ptype == "FeaturePlot") {
              my_guides <- guides(
                color = guide_colorbar(direction = leg_dir, title.position = "top", order = 1),
                fill = guide_colorbar(direction = leg_dir, title.position = "top", order = 1)
              )
            } else if (ptype == "ViolinPlot") {
              my_guides <- guides(fill = guide_legend(direction = leg_dir, override.aes = list(size = 4), title.position = "top"))
            } else if (ptype == "DimPlot") {
              my_guides <- guides(color = guide_legend(direction = leg_dir, override.aes = list(size = 4), title.position = "top"))
            } else if (ptype == "DotPlot") {
              my_guides <- guides(
                color = guide_colorbar(direction = leg_dir, title.position = "top", order = 1),
                size = guide_legend(direction = leg_dir, title.position = "top", order = 2)
              )
            }
            
            # Apply guides and theme using correct operator based on object class
            th <- theme(legend.position = legend_pos)
            
            if (inherits(p, "patchwork")) {
              # For patchwork, use & to apply to all parts and ensure guides are collected
              if (!is.null(my_guides)) p <- p & my_guides
              p <- p & th
              p <- p + patchwork::plot_layout(guides = "collect")
            } else {
              # For standard ggplot, use +
              if (!is.null(my_guides)) p <- p + my_guides
              p <- p + th
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
        colourpicker::colourInput(ns(paste0("col_", l)), paste("Color", l), value="gray")
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
      colourpicker::colourInput(paste0("hm_col_", l), paste("Color", l), value = "#808080")
    })
  })
  

  
  # File Upload Logic


  
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
  # --- NMF Analysis ---
  nmf_results <- reactiveValues(model = NULL, feature_loadings = NULL)

  observeEvent(input$run_nmf, {
    req(seurat_obj())
    obj <- seurat_obj()
    
    # Check dependencies
    if (!requireNamespace("NMF", quietly = TRUE)) {
      showNotification("NMF package is not installed.", type="error")
      return()
    }
    
    withProgress(message = 'Running NMF...', detail = 'Preparing data', value = 0.1, {
      
      # NMF requires non-negative data. Use counts (not data which may be log-normalized)
      # Check if we have variable features - if not, compute them automatically
      if (length(VariableFeatures(obj)) == 0) {
        showNotification("No variable features found. Computing automatically...", type = "message", duration = 3)
        
        # Clean data layer first by replacing NaN/Inf values
        data_layer <- GetAssayData(obj, layer = "data")
        data_layer[is.na(data_layer) | is.infinite(data_layer)] <- 0
        obj <- SetAssayData(obj, layer = "data", new.data = data_layer)
        
        # Now FindVariableFeatures should work
        obj <- tryCatch({
          FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
        }, error = function(e) {
          showNotification(paste("FindVariableFeatures error:", e$message), type = "warning", duration = 5)
          obj
        })
        
        # Check if we got variable features
        if (length(VariableFeatures(obj)) == 0) {
          showNotification("Using top 2000 expressed genes...", type = "warning", duration = 3)
          gene_means <- Matrix::rowMeans(data_layer)
          gene_means[is.na(gene_means) | is.infinite(gene_means)] <- 0
          top_genes <- names(sort(gene_means, decreasing = TRUE)[1:min(2000, sum(gene_means > 0))])
          VariableFeatures(obj) <- top_genes
        }
        
        seurat_obj(obj)
      }
      
      # Get RAW counts for NMF (not normalized data)
      # H5AD conversion may use different layer names
      mat <- NULL
      
      # Try to get raw counts from various possible locations
      for (layer_name in c("counts.raw", "raw.counts", "counts", "data")) {
        temp <- tryCatch({
          GetAssayData(obj, layer = layer_name)
        }, error = function(e) NULL)
        
        if (!is.null(temp) && all(temp >= 0, na.rm = TRUE)) {
          message("Using layer: ", layer_name)
          mat <- temp
          break
        }
      }
      
      # If still NULL, try accessing raw data directly
      if (is.null(mat)) {
        mat <- tryCatch({
          obj@assays$RNA@layers$counts
        }, error = function(e) {
          GetAssayData(obj, layer = "data")
        })
      }
      
      mat <- mat[VariableFeatures(obj), ]
      
      # Clean the matrix: replace NaN/Inf with 0, and ensure non-negative
      mat[is.na(mat) | is.infinite(mat)] <- 0
      
      # Check if values are extremely small (normalized/scaled data)
      # NMF has issues with very small decimal values
      max_val <- max(mat, na.rm = TRUE)
      if (max_val > 0 && max_val < 1) {
        showNotification("Data appears normalized. Scaling up for NMF...", type = "message", duration = 3)
        # Scale up to reasonable range (multiply by 1000)
        mat <- mat * 1000
      }
      
      # H5AD files sometimes have negative values even in counts - force non-negative
      if (any(mat < 0, na.rm = TRUE)) {
        showNotification("Converting negative values to positive (abs).", type = "warning", duration = 3)
        mat <- abs(mat)
      }
      
      # Remove zero rows
      mat <- mat[rowSums(mat) > 0, ]
      
      # Convert to dense matrix early (sparse matrices can cause issues with NMF)
      if (inherits(mat, "dgCMatrix") || inherits(mat, "sparseMatrix")) {
        showNotification("Converting sparse matrix to dense for NMF...", type = "message", duration = 2)
        mat <- as.matrix(mat)
      }
      
      # Validate matrix after filtering
      if (nrow(mat) == 0) {
        showNotification("All genes have zero expression after filtering. Cannot run NMF.", type = "error")
        return()
      }
      
      if (nrow(mat) < input$nmf_k) {
        showNotification(paste("Not enough genes (", nrow(mat), ") for", input$nmf_k, "factors. Please reduce k."), type = "error")
        return()
      }
      
      incProgress(0.2, detail = paste("Processing", nrow(mat), "genes x", ncol(mat), "cells"))
      
      # NMF package requires dense matrix
      mat <- as.matrix(mat)
      
      # Scale data up if values are too small (prevents numerical issues)
      # Multiply by 1000 and round to avoid floating point precision issues
      mat <- round(mat * 1000)
      
      # Ensure strictly non-negative after rounding
      mat[mat < 0] <- 0
      
      # Final validation before NMF
      if (any(is.na(mat))) {
        showNotification("Matrix contains NA values. Cleaning...", type = "warning", duration = 2)
        mat[is.na(mat)] <- 0
      }
      if (any(is.infinite(mat))) {
        showNotification("Matrix contains Inf values. Cleaning...", type = "warning", duration = 2)
        mat[is.infinite(mat)] <- 0
      }
      if (any(mat < 0)) {
        showNotification("Matrix contains negative values. Converting to absolute...", type = "warning", duration = 2)
        mat <- abs(mat)
      }
      
      incProgress(0.3, detail = paste("Computing", input$nmf_k, "factors (10-20 min)..."))
      
      k <- input$nmf_k
      
      # Log matrix details for debugging
      message("=== NMF Debug Info ===")
      message("Matrix class: ", paste(class(mat), collapse = ", "))
      message("Matrix type: ", typeof(mat))
      message("Is matrix: ", is.matrix(mat))
      message("Dimensions: ", nrow(mat), " x ", ncol(mat))
      message("Value range: ", min(mat, na.rm = TRUE), " to ", max(mat, na.rm = TRUE))
      message("k value: ", k)
      message("method: snmf/l (fast)")
      message("======================")
      
      # Run NMF
      res <- tryCatch({
        NMF::nmf(mat, rank = k, method = "snmf/r", seed = 123456)
      }, error = function(e) {
        showNotification(paste("NMF error:", e$message), type = "error", duration = 10)
        message("NMF failed with error: ", e$message)
        message("Matrix dimensions: ", nrow(mat), " x ", ncol(mat))
        message("Matrix range: ", min(mat, na.rm = TRUE), " to ", max(mat, na.rm = TRUE))
        message("k value: ", k)
        return(NULL)
      })
      
      incProgress(0.8, detail="Processing results...")
      
      # Check if NMF succeeded
      if (is.null(res)) {
        showNotification("NMF computation failed. Check console for details.", type = "error", duration = 10)
        return()
      }
      
      W <- NMF::basis(res)
      H <- NMF::coef(res)
      
      factor_names <- paste0("Factor_", 1:k)
      colnames(W) <- factor_names
      rownames(H) <- factor_names
      colnames(H) <- colnames(mat)
      
      nmf_results$model <- res
      nmf_results$feature_loadings <- W
      
      embeddings <- t(H)
      rownames(embeddings) <- colnames(mat)
      colnames(embeddings) <- factor_names
      
      obj[["nmf"]] <- CreateDimReducObject(embeddings = embeddings, loadings = W, key = "NMF_", assay = DefaultAssay(obj))
      
      seurat_obj(obj)
      showNotification("NMF Analysis Complete!", type="message")
      
      updateSelectInput(session, "nmf_factor", choices = factor_names, selected = factor_names[1])
    })
  })
  
  output$nmf_done <- reactive({
    return(!is.null(nmf_results$model))
  })
  outputOptions(output, "nmf_done", suspendWhenHidden = FALSE)
  
  output$nmf_heatmap <- renderPlot({
    req(nmf_results$feature_loadings)
    W <- nmf_results$feature_loadings
    top_genes <- unique(as.vector(apply(W, 2, function(x) names(sort(x, decreasing=TRUE))[1:10])))
    W_sub <- W[top_genes, ]
    df <- reshape2::melt(W_sub)
    colnames(df) <- c("Gene", "Factor", "Weight")
    
    ggplot(df, aes(x=Factor, y=Gene, fill=Weight)) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=45, hjust=1), axis.text.y = element_text(size=8)) +
      labs(title = "Gene Programs (Basis Matrix)")
  })
  
  output$nmf_feature_plot <- renderPlot({
    req(seurat_obj(), input$nmf_factor)
    obj <- seurat_obj()
    
    # Check if NMF has been run
    if (!"nmf" %in% names(obj@reductions)) {
      plot.new()
      text(0.5, 0.5, "Run NMF first to see factor scores.", cex = 1.5)
      return()
    }
    
    nmf_data <- Embeddings(obj, "nmf")
    factor_col <- input$nmf_factor
    
    # Validate factor exists
    if (is.null(factor_col) || !factor_col %in% colnames(nmf_data)) {
      plot.new()
      text(0.5, 0.5, "Select a factor from the sidebar.", cex = 1.5)
      return()
    }
    
    obj$current_factor <- nmf_data[, factor_col]
    obj$current_factor <- nmf_data[, factor_col]
    
    red <- input$nmf_reduction
    if (is.null(red) || !red %in% Reductions(obj)) {
      avail <- Reductions(obj)
      if (length(avail) > 0) red <- avail[1] else red <- NULL
    }
    
    if (is.null(red)) {
       plot.new(); text(0.5, 0.5, "No reductions (UMAP/tSNE) found."); return()
    }
    
    FeaturePlot(obj, features = "current_factor", reduction = red, pt.size = input$nmf_pt_size) +
      scale_color_viridis_c(option = tolower(input$nmf_color)) +
      labs(title = paste(factor_col, "Scores"))
  })
  
  output$nmf_gene_table <- DT::renderDT({
    req(nmf_results$feature_loadings, input$nmf_factor)
    W <- nmf_results$feature_loadings
    col <- input$nmf_factor
    scores <- W[, col]
    sorted_scores <- sort(scores, decreasing=TRUE)
    
    # Show top 50 genes for the selected factor
    top_n <- min(50, length(sorted_scores))
    df <- data.frame(
      Rank = 1:top_n,
      Gene = names(sorted_scores)[1:top_n], 
      Weight = round(sorted_scores[1:top_n], 4),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(df, 
      options = list(
        pageLength = 20,
        dom = 'ftip',
        ordering = TRUE
      ),
      rownames = FALSE,
      caption = paste("Top genes defining", col)
    )
  })
  
  # New: Multi-factor FeaturePlot
  output$nmf_all_factors_plot <- renderPlot({
    req(seurat_obj())
    obj <- seurat_obj()
    if (!"nmf" %in% names(obj@reductions)) {
      return(NULL)
    }
    
    nmf_data <- Embeddings(obj, "nmf")
    factors <- colnames(nmf_data)
    
    # Add all factors to metadata for plotting
    for (f in factors) {
      obj@meta.data[[f]] <- nmf_data[, f]
    }
    
    # Determine reduction
    red <- input$nmf_reduction
    if (!red %in% Reductions(obj)) {
      # Fallback if selected reduction doesn't exist
      avail <- Reductions(obj)
      if (length(avail) > 0) red <- avail[1] else return(NULL)
    }
    
    # Create multi-panel FeaturePlot
    ncols <- ceiling(sqrt(length(factors)))
    
    # Use SCpubr if available, otherwise standard Seurat
    if (has_scpubr) {
      SCpubr::do_FeaturePlot(obj, features = factors, 
                            reduction = red,
                            pt.size = input$nmf_pt_size,
                            ncol = ncols,
                            legend.position = "bottom")
    } else {
      FeaturePlot(obj, features = factors, 
                 reduction = red,
                 pt.size = input$nmf_pt_size,
                 ncol = ncols) &
        scale_color_viridis_c(option = tolower(input$nmf_color))
    }
  })
  
  # Export handlers
  nmf_heatmap_plot <- reactive({
    req(nmf_results$feature_loadings)
    W <- nmf_results$feature_loadings
    top_genes <- unique(as.vector(apply(W, 2, function(x) names(sort(x, decreasing=TRUE))[1:10])))
    W_sub <- W[top_genes, ]
    df <- reshape2::melt(W_sub)
    colnames(df) <- c("Gene", "Factor", "Weight")
    
    ggplot(df, aes(x=Factor, y=Gene, fill=Weight)) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=45, hjust=1), axis.text.y = element_text(size=8)) +
      labs(title = "Gene Programs (Basis Matrix)")
  })
  
  output$download_nmf_heatmap <- downloadHandler(
    filename = function() {
      paste0("NMF_Heatmap_", Sys.Date(), ".", input$nmf_export_format)
    },
    content = function(file) {
      req(nmf_heatmap_plot())
      ggsave(file, nmf_heatmap_plot(),
             width = input$nmf_export_width,
             height = input$nmf_export_height,
             device = input$nmf_export_format)
    }
  )
  
  nmf_featureplot_reactive <- reactive({
    req(seurat_obj())
    obj <- seurat_obj()
    if (!"nmf" %in% names(obj@reductions)) return(NULL)
    
    nmf_data <- Embeddings(obj, "nmf")
    factors <- colnames(nmf_data)
    
    for (f in factors) {
      obj@meta.data[[f]] <- nmf_data[, f]
    }
    
    red <- if ("umap" %in% Reductions(obj)) "umap" else "pca"
    ncols <- ceiling(sqrt(length(factors)))
    
    if (has_scpubr) {
      SCpubr::do_FeaturePlot(obj, features = factors, 
                            reduction = red,
                            pt.size = input$nmf_pt_size,
                            ncol = ncols,
                            legend.position = "bottom")
    } else {
      FeaturePlot(obj, features = factors, 
                 reduction = red,
                 pt.size = input$nmf_pt_size,
                 ncol = ncols) &
        scale_color_viridis_c(option = tolower(input$nmf_color))
    }
  })
  
  # --- Multimodal Logic ---
  
  # Reactive for Co-expression Plot to support download
  adt_coexp_plot <- reactive({
    req(seurat_obj(), input$adt_feature_x, input$rna_feature_y)
    obj <- seurat_obj()
    
    # Get Data
    df <- FetchData(obj, vars = c(input$adt_feature_x, input$rna_feature_y))
    colnames(df) <- c("Protein", "RNA")
    
    # Theme settings
    font_size <- if(!is.null(input$adt_font_size)) input$adt_font_size else 12
    pt_size <- if(!is.null(input$adt_pt_size)) input$adt_pt_size else 1
    
    p <- ggplot(df, aes(x=Protein, y=RNA)) +
      geom_point(alpha=0.6, size=pt_size, color="#2c3e50") +
      stat_density_2d(color="#e74c3c") +
      theme_minimal(base_size = font_size) +
      labs(title = paste(input$adt_feature_x, "vs", input$rna_feature_y),
           x = input$adt_feature_x, y = input$rna_feature_y) +
      geom_smooth(method="lm", color="blue", se=FALSE)
      
    return(p)
  })

  output$adt_rna_coexpression <- renderPlot({
    req(adt_coexp_plot())
    adt_coexp_plot()
  })
  
  # Reactive for FeaturePlot (ADT)
  output$adt_feature_plot <- renderPlot({
    req(seurat_obj(), input$adt_feature)
    obj <- seurat_obj()
    
    # Palettes
    pal <- input$adt_palette
    if (is.null(pal)) pal <- "viridis"
    
    # Map to SCpubr/Viridis options if needed, or use base logic
    # Simplified approach: Use FeaturePlot and apply scale
    p <- FeaturePlot(obj, features = input$adt_feature, min.cutoff = "q9") +
      theme_void() +
      labs(title = input$adt_feature)
      
    # Apply palette
    if (pal %in% c("viridis", "magma", "plasma", "inferno", "cividis")) {
      p <- p + scale_color_viridis_c(option = pal)
    } else if (pal %in% c("RdBu", "RdYlBu", "Spectral")) {
      p <- p + scale_color_distiller(palette = pal)
    }
    
    # Apply sizing
    pt_size <- if(!is.null(input$adt_pt_size)) input$adt_pt_size else 1
    # Seurat FeaturePlot uses 'pt.size' argument, but since we already generated p...
    # We can try to update the layer or just regenerate with pt.size if possible. 
    # Actually FeaturePlot takes pt.size. Let's regenerate.
    
    p <- FeaturePlot(obj, features = input$adt_feature, min.cutoff = "q9", pt.size = pt_size) +
         theme_void(base_size = if(!is.null(input$adt_font_size)) input$adt_font_size else 12) +
         labs(title = input$adt_feature)

    if (pal %in% c("viridis", "magma", "plasma", "inferno", "cividis")) {
      p <- p + scale_color_viridis_c(option = pal)
    } else {
        p <- p + scale_color_distiller(palette = pal)
    }
    p
  })
  
  # Download Handlers for ADT
  output$download_adt_coexp <- downloadHandler(
    filename = function() { paste0("Coexpression_", input$adt_feature_x, "_vs_", input$rna_feature_y, ".", input$adt_export_format) },
    content = function(file) {
      req(adt_coexp_plot())
      ggsave(file, adt_coexp_plot(), width = input$adt_export_width, height = input$adt_export_height)
    }
  )

  output$download_nmf_featureplot <- downloadHandler(
    filename = function() {
      paste0("NMF_FeaturePlot_AllFactors_", Sys.Date(), ".", input$nmf_export_format)
    },
    content = function(file) {
      req(nmf_featureplot_reactive())
      ggsave(file, nmf_featureplot_reactive(),
             width = input$nmf_export_width,
             height = input$nmf_export_height,
             device = input$nmf_export_format)
    }
  )
  
  output$download_nmf_table <- downloadHandler(
    filename = function() {
      paste0("NMF_TopGenes_", input$nmf_factor, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(nmf_results$feature_loadings, input$nmf_factor)
      W <- nmf_results$feature_loadings
      col <- input$nmf_factor
      scores <- W[, col]
      req(seurat_obj(), input$nmf_gene_table_rows_all)
      # Get the correct indices from DT (filtered/sorted)
      idx <- input$nmf_gene_table_rows_all
      
      # Reconstruct the table logic (simplified)
      res <- nmf_results()
      req(res)
      W <- basis(res)
      factors <- input$nmf_factor
      if (is.null(factors)) factors <- colnames(W)[1]
      
      scores <- W[, factors]
      df <- data.frame(Gene = names(sort(scores, decreasing=TRUE)), 
                      Weight = sort(scores, decreasing=TRUE))
      write.csv(df, file, row.names = FALSE)
    }
  )

  # --- Gating Logic ---
  
  output$gating_plot <- renderPlotly({
    req(seurat_obj(), input$gate_feature_x, input$gate_feature_y)
    
    obj <- seurat_obj()
    feat_x <- input$gate_feature_x
    feat_y <- input$gate_feature_y
    
    # Fetch data (check ADT first then RNA)
    # Helper to get data
    get_data <- function(f, o) {
      if (f %in% rownames(o)) {
        FetchData(o, vars = f)[,1]
      } else {
        # Try finding in ADT
        if ("ADT" %in% names(o@assays) && f %in% rownames(o[["ADT"]])) {
           FetchData(o, vars = paste0("ADT_", f))[,1] # Seurat often prefixes
           # If fails, try raw
           d <- GetAssayData(o, assay="ADT", layer="data")[f, ]
           return(d)
        }
        return(NULL)
      }
    }
    
    # Simple FetchData is safer
    df <- FetchData(obj, vars = c(feat_x, feat_y))
    colnames(df) <- c("X", "Y")
    df$CellID <- rownames(df)
    
    plot_ly(df, x = ~X, y = ~Y, key = ~CellID, type = "scatter", mode = "markers",
            marker = list(size = 5, opacity = 0.6)) %>%
      layout(
        title = paste(feat_x, "vs", feat_y),
        xaxis = list(title = feat_x),
        yaxis = list(title = feat_y),
        dragmode = "lasso"
      ) %>%
      event_register("plotly_selected")
  })
  
  output$gating_info <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) "Select cells on the plot using Lasso or Box tool."
    else paste("Selected", nrow(d), "cells.")
  })
  
  observeEvent(input$gate_cells, {
    req(seurat_obj())
    d <- event_data("plotly_selected")
    if (is.null(d)) {
      showNotification("No cells selected!", type = "warning")
      return()
    }
    
    selected_cells <- d$key
    if (length(selected_cells) > 0) {
      obj_sub <- subset(seurat_obj(), cells = selected_cells)
      seurat_obj(obj_sub)
      showNotification(paste("Gated to", length(selected_cells), "cells."), type = "message")
    }
  })
  
  observeEvent(input$reset_gating, {
    req(original_obj())
    seurat_obj(original_obj())
    showNotification("Reset to original dataset.", type = "message")
  })

} # End Server
shinyApp(ui, server)
