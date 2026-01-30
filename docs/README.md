# Documentation

This directory contains technical documentation for the Seurat Interactive Visualizer project.

## Contents

### Development Documentation

- **[REFACTORING_WALKTHROUGH.md](REFACTORING_WALKTHROUGH.md)** - Detailed walkthrough of the modular refactoring work, including bug fixes, code metrics, and testing instructions
- **[REFACTORING_PLAN.md](REFACTORING_PLAN.md)** - Implementation plan for reducing code complexity through modular architecture

### Setup & Deployment

- **[RENV_SETUP.md](../RENV_SETUP.md)** - Instructions for setting up renv for reproducible package management
- **[DOCKER_INSTRUCTIONS.md](../DOCKER_INSTRUCTIONS.md)** - Docker build and deployment instructions

### Code Organization

The refactoring introduced modular utility files:

- `color_utils.R` - Unified color palette system (viridis, MetBrewer, RColorBrewer)
- `plot_cluster_distribution.R` - Self-contained ClusterDistrBar plotting module

Future modules (planned):
- `plot_dimension_reduction.R` - DimPlot module
- `plot_feature.R` - FeaturePlot module
- `plot_violin.R` - ViolinPlot module
- `plot_dot.R` - DotPlot module

## Contributing

When making significant changes:

1. Update or create documentation in this directory
2. Follow the modular pattern established in the refactoring
3. Use renv to track package dependencies
4. Test changes in Docker before committing
5. Write clear commit messages following conventional commits format

## Version History

- **2026-01-29** - Initial refactoring: unified color system and ClusterDistrBar module
  - Commit: `75b79f4` - feat: Add unified color utilities module
  - Fixed ClusterDistrBar manual color picker bug
  - Reduced code complexity by ~70%
