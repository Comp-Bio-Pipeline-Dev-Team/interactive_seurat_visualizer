# Documentation

This directory is the **single source of truth** for all project documentation. All planning artifacts, walkthroughs, setup guides, and session prompts are kept here and version-controlled alongside the code.

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
├── # Feature Planning (active & completed)
├── SIMPLIFICATION_PLAN_03292026.md    ← Plan: remove plotly/DE/ORA tabs (in progress)
├── LANDING_PAGE_IMPLEMENTATION_PLAN.md← Plan: landing page with metadata mgmt
├── LANDING_PAGE_STATUS.md             ← Status: landing page feature progress
├── REFACTORING_PLAN.md                ← Plan: modular architecture refactoring
├── MODULARIZATION_PLAN_V2.md          ← Plan: V2 modularization
│
├── # Walkthroughs (completed work)
├── REFACTORING_WALKTHROUGH.md         ← Walkthrough: modular refactoring results
├── MODULARIZATION_WALKTHROUGH_V2.md   ← Walkthrough: V2 modularization
├── MODULARIZATION_WALKTHROUGH_FINAL.md← Walkthrough: Final modularization state
│
├── # Testing & Reference
├── TESTING_LANDING_PAGE.md            ← Testing guide: landing page features
├── PACKAGE_AUDIT.md                   ← R package dependency audit
├── RENV_LOCK_GENERATION.md            ← How renv.lock was generated via Docker
│
└── antigravity_prompting/             ← Full AI session transcripts (auto-saved)
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

> **All planning artifacts, implementation walkthroughs, and session prompts go in `docs/`. Always commit them with (or before) the code changes they describe.**

| Doc type | Naming convention | When to create |
|---|---|---|
| Implementation plan | `*_PLAN_MMDDYYYY.md` | Before starting significant changes |
| Walkthrough | `*_WALKTHROUGH_MMDDYYYY.md` | After completing changes |
| AI session prompt | `docs/antigravity_prompting/Topic_antigravity_prompting_MMDDYYYY.md` | Auto-saved by Antigravity |

Workflow:
1. **Plan** → create `docs/*_PLAN_MMDDYYYY.md`, commit it
2. **Execute** → make code changes
3. **Document** → create `docs/*_WALKTHROUGH_MMDDYYYY.md`, commit with code
4. **Update** → this `README.md` feature history table

---

## 📋 Feature History

| Date | Feature | Status | Plan | Walkthrough |
|------|---------|--------|------|-------------|
| 2026-01-28 | Modular architecture refactoring | ✅ Done | [REFACTORING_PLAN.md](REFACTORING_PLAN.md) | [REFACTORING_WALKTHROUGH.md](REFACTORING_WALKTHROUGH.md) |
| 2026-01-29 | V2 modularization (plot modules) | ✅ Done | [MODULARIZATION_PLAN_V2.md](MODULARIZATION_PLAN_V2.md) | [MODULARIZATION_WALKTHROUGH_FINAL.md](MODULARIZATION_WALKTHROUGH_FINAL.md) |
| 2026-01-30 | Landing page & metadata management | ✅ Done | [LANDING_PAGE_IMPLEMENTATION_PLAN.md](LANDING_PAGE_IMPLEMENTATION_PLAN.md) | [LANDING_PAGE_STATUS.md](LANDING_PAGE_STATUS.md) |
| 2026-02-03 | UI enhancements (accordion, sidebar, legends, axes) | ✅ Done | — | [antigravity_prompting/](antigravity_prompting/) |
| 2026-03-29 | Simplification: remove plotly/DE/ORA | ⏳ In progress | [SIMPLIFICATION_PLAN_03292026.md](SIMPLIFICATION_PLAN_03292026.md) | — |

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
| `renv.lock` | **Pinned R package versions (source of truth for Docker)** |
| `Dockerfile` | Production Docker image (builds from renv.lock) |

