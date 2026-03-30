# Documentation

This directory is the **single source of truth** for all project documentation. All planning artifacts, walkthroughs, setup guides, and session logs are kept here and version-controlled alongside the code.

---

## 📁 Directory Structure

```
docs/
├── README.md                          ← This file (documentation index)
│
├── # Setup & Deployment
├── QUICKSTART.md                      ← How to build and run the app quickly
├── DOCKER_COMMANDS.md                 ← Docker build, run, and management commands
├── DOCKER_INSTRUCTIONS.md             ← Docker deployment details and volume mounting
├── RENV_SETUP.md                      ← renv package management and lockfile workflow
│
├── # Feature Implementation Docs
├── REFACTORING_PLAN.md                ← Plan: modular architecture refactoring
├── REFACTORING_WALKTHROUGH.md         ← Walkthrough: modular refactoring results
├── MODULARIZATION_PLAN_V2.md          ← Plan: V2 modularization
├── MODULARIZATION_WALKTHROUGH_V2.md   ← Walkthrough: V2 modularization
├── MODULARIZATION_WALKTHROUGH_FINAL.md← Walkthrough: Final modularization state
├── LANDING_PAGE_IMPLEMENTATION_PLAN.md← Plan: landing page with metadata mgmt
├── LANDING_PAGE_STATUS.md             ← Status: landing page feature progress
├── TESTING_LANDING_PAGE.md            ← Testing guide: landing page features
├── PACKAGE_AUDIT.md                   ← R package dependency audit
├── RENV_LOCK_GENERATION.md            ← How renv.lock was generated via Docker
│
└── antigravity_prompting/             ← AI-assisted session logs
    ├── Refactoring_Seurat_App_antigravity_prompting_01302026.md
    └── Implementing_UI_Enhancements_antigravity_prompting_02032026.md
```

---

## 🐳 Docker & Reproducibility

The Docker image is built **from `renv.lock`** to ensure fully reproducible package environments:

```bash
# Build image (uses renv.lock to restore all packages)
sudo docker build -t seurat-shiny-app .

# Run container
sudo docker run -d -p 3838:3838 --name seurat-app seurat-shiny-app

# Access app
open http://localhost:3838
```

### Dependency Chain (version controlled)
```
renv.lock  ──►  Dockerfile  ──►  Docker image  ──►  Running container
   ↑
renv::snapshot() to update
```

**When adding new R packages:**
1. Add the package to your R code
2. Run `renv::snapshot()` inside the container to update `renv.lock`
3. Rebuild the Docker image
4. Commit both the code change and the updated `renv.lock`

See [`RENV_SETUP.md`](RENV_SETUP.md) and [`DOCKER_COMMANDS.md`](DOCKER_COMMANDS.md) for full details.

---

## 🗂️ Documentation Convention

> **All planning artifacts, implementation walkthroughs, and session notes go in `docs/`.**

When making significant changes:
1. Create or update a `*_PLAN.md` in `docs/` before starting
2. Create or update a `*_WALKTHROUGH.md` in `docs/` after completing
3. Session AI-prompting logs go in `docs/antigravity_prompting/`
4. Commit docs changes **with** the code changes they describe

---

## 📋 Feature History

| Date | Feature | Plan | Walkthrough |
|------|---------|------|-------------|
| 2026-01-28 | Modular architecture refactoring | [REFACTORING_PLAN.md](REFACTORING_PLAN.md) | [REFACTORING_WALKTHROUGH.md](REFACTORING_WALKTHROUGH.md) |
| 2026-01-29 | V2 modularization (plot modules) | [MODULARIZATION_PLAN_V2.md](MODULARIZATION_PLAN_V2.md) | [MODULARIZATION_WALKTHROUGH_FINAL.md](MODULARIZATION_WALKTHROUGH_FINAL.md) |
| 2026-01-30 | Landing page & metadata management | [LANDING_PAGE_IMPLEMENTATION_PLAN.md](LANDING_PAGE_IMPLEMENTATION_PLAN.md) | [LANDING_PAGE_STATUS.md](LANDING_PAGE_STATUS.md) |
| 2026-02-03 | UI enhancements (accordion, sidebar, legends, axes) | — | See session log |

---

## 🏗️ Code Organization

| File | Purpose |
|------|---------|
| `app.R` | Main Shiny application |
| `color_utils.R` | Unified color palette system |
| `plot_dimension_reduction.R` | DimPlot module |
| `plot_feature.R` | FeaturePlot module |
| `plot_violin.R` | ViolinPlot module |
| `plot_dot.R` | DotPlot module |
| `plot_cluster_distribution.R` | ClusterDistrBar module |
| `ui_landing_page.R` | Landing page UI |
| `server_landing_page.R` | Landing page server logic |
| `setup.R` | App initialization helpers |
| `enrichment_backend.R` | Gene enrichment analysis backend |
| `renv.lock` | **Pinned R package versions (source of truth for Docker)** |
| `Dockerfile` | Production Docker image (builds from renv.lock) |
