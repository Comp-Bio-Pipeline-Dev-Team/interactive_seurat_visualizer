# Simplification: Remove Plotly, DE Tab, and Pathway Enrichment Tab

Remove three features to reduce app complexity, dependencies, and Docker build time.

## What's Being Removed

| Feature | Where It Lives | Scope |
|---|---|---|
| Interactive plotly toggle | `app.R` UI + server | Remove `plot_style` selector + `renderPlotly` in volcano; swap `plotlyOutput` → `plotOutput` for volcano |
| Differential Expression tab | `app.R` UI (lines 312–358) + server (lines 680–807) | Remove entire tab + all server observers/reactives |
| Pathway Enrichment tab | `app.R` UI (lines 361–445) + server (lines 808–1264) + `enrichment_backend.R` | Remove tab, all server logic, and the entire backend file |

---

## Proposed Changes

### `app.R`

#### Libraries (lines 1–26)
Remove `library(plotly)` and all optional enrichment libraries:
```diff
- library(plotly)
- has_enrichment <- requireNamespace("clusterProfiler"...) ...
- if (has_enrichment) { library(clusterProfiler); library(enrichplot); ... }
```

#### `plotControlUI` function (line 95)
Remove the `"Interactive (plotly)"` option from the Plot Style selector:
```diff
- selectInput(ns("plot_style"), "Plot Style", choices = c("Static (ggplot2)", "Interactive (plotly)"))
+ # Remove Plot Style selector entirely (all plots are now static ggplot2)
```

#### UI: Differential Expression tab (lines 312–358)
Delete entire `tabPanel("Differential Expression", ...)` block.

#### UI: Pathway Enrichment tab (lines 361–445)
Delete entire `tabPanel("Pathway Enrichment", ...)` block.

#### Server: DE reactive values & logic (lines 680–807)
Delete:
- `de_results <- reactiveVal(NULL)` 
- All `observe/observeEvent` blocks for `de_group_by`, `run_de`
- `output$de_table`, `volcano_plot`, `output$de_volcano`, `output$download_volcano`, `output$download_de`

#### Server: Pathway enrichment logic (lines 808–1264)
Delete entire block between `# ===== PATHWAY ENRICHMENT LOGIC =====` and `# ===== END PATHWAY ENRICHMENT LOGIC =====`

#### Server: `generate_plot` function
Remove `plot_style` branch — all plots are ggplot2 objects returned directly. The `output[[ns("plot_style_ui")]]` renderUI already renders only `"Standard"`. Remove that `renderUI` too.

---

### `enrichment_backend.R`

Delete this file entirely (no longer used).

---

### `renv.lock` and `Dockerfile`

Packages to remove from `renv.lock` (no longer needed):
- `plotly` + its dependencies (`crosstalk`, `htmlwidgets`, `lazyeval`, `digest` if plotly-only)
- `clusterProfiler`, `enrichplot`, `fgsea`, `msigdbr`, `DOSE`, `igraph`
- `org.Hs.eg.db`, `org.Mm.eg.db`, `AnnotationDbi` (only used in enrichment + Ensembl conversion)

> [!IMPORTANT]
> `AnnotationDbi` / `org.*.eg.db` are also used by the **Ensembl ID conversion** feature on the landing page. These must be **kept** if we keep that feature, or removed only if we remove/stub that feature too.

> [!NOTE]  
> `renv.lock` is regenerated inside Docker via `renv::snapshot()`. The workflow is: remove library calls → rebuild image → snapshot updated lockfile → commit. The current Dockerfile already runs `renv::restore()` from `renv.lock`, so we need to first update the R code, then regenerate the lockfile.

---

## Open Questions

> [!IMPORTANT]
> **Ensembl ID conversion** (on the landing page) uses `AnnotationDbi`, `org.Hs.eg.db`, and `org.Mm.eg.db`. These overlap with the enrichment packages being removed. Should we:
> - **Keep** those 3 packages (they're relatively lean, and the conversion feature is useful)
> - **Remove** the conversion feature as well to further simplify?

---

## Verification Plan

### Automated / In-Docker
1. Rebuild Docker image after code changes
2. Run container, upload `test_seurat.rds`
3. Verify:
   - Landing page works, Proceed to App works
   - Visualization tab: all 5 plot types (DimPlot, FeaturePlot, ViolinPlot, DotPlot, ClusterDistrBar) render correctly
   - Heatmap tab works
   - DE and Pathway Enrichment tabs are **gone** from navbar
   - No plotly-related errors in console

### Manual
- Confirm sidebar plot controls show no "Interactive (plotly)" option
- Confirm export/download still works for ggplot2 plots
