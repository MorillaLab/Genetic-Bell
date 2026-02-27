# Genetic-Bell R Package Dependencies
# Run this script in R to install all required packages:
#   source("requirements.R")

required_packages <- c(
  # Network analysis
  "igraph",         # Gene interaction network construction and community detection

  # Statistical / mathematical methods
  "MASS",           # MDS (multidimensional scaling), multivariate analysis
  "cluster",        # k-means clustering, silhouette analysis
  "mixtools",       # EM algorithm for mixture models (peak detection)
  "stats",          # Covariance matrix, PCA (base R)

  # GWAS / genomics data handling
  "data.table",     # Fast reading of large GWAS matrices
  "Matrix",         # Sparse matrix operations

  # Visualisation
  "ggplot2",        # Publication-quality figures
  "corrplot",       # Covariance / correlation matrix heatmaps
  "ggraph",         # Network visualisation (ggplot2 extension)
  "plotly",         # Interactive 3D CHSH space plots
  "rgl",            # 3D polyhedra visualisation (CHSH space)
  "patchwork",      # Multi-panel figure composition

  # Utilities
  "tidyverse",      # Data manipulation (dplyr, tidyr, readr)
  "scales",         # Axis scaling helpers
  "RColorBrewer"    # Colour palettes for figures
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  install.packages(new_packages, repos = "https://cloud.r-project.org")
  cat("Installed:", paste(new_packages, collapse = ", "), "\n")
} else {
  cat("All required packages are already installed.\n")
}
