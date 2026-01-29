# Pathway Enrichment Backend Logic
# This file contains the server-side logic for enrichment analysis
# To be inserted into app.R after the download_de handler

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
      genes <- df$gene
      gene_fc <- setNames(df$avg_log2FC, df$gene)
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
        enrichment_results(result)
        showNotification(paste("Found", nrow(result@result), "enriched pathways"), 
                        type = "message")
      } else {
        showNotification("No significant enrichment found", type = "warning")
        enrichment_results(NULL)
      }
    })
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
  
  # Dot plot
  output$enrichment_dotplot <- renderPlot({
    req(enrichment_results())
    tryCatch({
      dotplot(enrichment_results(), showCategory = 20) + 
        theme(axis.text.y = element_text(size = 10))
    }, error = function(e) {
      ggplot() + annotate("text", x = 0, y = 0, 
                         label = paste("Plot error:", e$message), 
                         size = 4, color = "red") + theme_void()
    })
  })
  
  # Bar plot
  output$enrichment_barplot <- renderPlot({
    req(enrichment_results())
    tryCatch({
      barplot(enrichment_results(), showCategory = 20) +
        theme(axis.text.y = element_text(size = 10))
    }, error = function(e) {
      ggplot() + annotate("text", x = 0, y = 0, 
                         label = paste("Plot error:", e$message), 
                         size = 4, color = "red") + theme_void()
    })
  })
  
  # Network plot
  output$enrichment_network <- renderPlot({
    req(enrichment_results())
    tryCatch({
      # Need to calculate pairwise similarities
      result_sim <- pairwise_termsim(enrichment_results())
      emapplot(result_sim, showCategory = 30)
    }, error = function(e) {
      ggplot() + annotate("text", x = 0, y = 0, 
                         label = paste("Network plot error:", e$message), 
                         size = 4, color = "red") + theme_void()
    })
  })
  
  # GSEA enrichment plots
  output$enrichment_gsea <- renderPlot({
    req(enrichment_results(), input$enrich_type == "gsea")
    tryCatch({
      # Show top 5 pathways
      n_pathways <- min(5, nrow(enrichment_results()))
      gseaplot2(enrichment_results(), geneSetID = 1:n_pathways, 
               pvalue_table = TRUE)
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
  
  # ===== END PATHWAY ENRICHMENT LOGIC =====
