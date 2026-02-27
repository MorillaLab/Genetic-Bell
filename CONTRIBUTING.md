# Contributing to Genetic-Bell

Thank you for your interest! This is a mathematically deep project — contributions that extend the CHSH space model, apply it to new gene pairs or diseases, or improve the R implementation are very welcome.

## 🐛 Reporting Bugs

Open a [GitHub Issue](https://github.com/MorillaLab/Genetic-Bell/issues) with:
- The R script or function where the error occurs
- Your R version (`R.version`) and package versions
- A minimal reproducible example
- The full error message / traceback

## 💡 Suggesting Features or Extensions

Open an issue tagged `enhancement`. Good examples:
- New CHSH inequality formulations (e.g., higher-order interactions, 3+ genes)
- Integration with eQTL data, Hi-C chromatin contact data
- Application to a new complex disease or SNP dataset
- Bayesian extension of the module permissivity model
- Improved 3D CHSH space visualisation

## 🔧 Submitting Code

1. Fork the repository and create a branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```
2. Install dependencies:
   ```r
   source("requirements.R")
   ```
3. Write clean, well-commented R code. For any new mathematical object, include:
   - A comment with the formal definition
   - A reference to the relevant paper (Bell 1964, CHSH 1969, or your extension)
4. Add a short example or test in the `model/` folder
5. Open a pull request against `main` with a clear mathematical and biological motivation

## 📐 Mathematical Conventions

- CHSH space coordinates follow Clauser-Horne-Shimony-Holt (1969)
- Tensor notation should match the paper (Kaba *et al.*, 2023)
- New module definitions should reference a GO term or pathway database (KEGG, Reactome)
- Variable names should be descriptive: `chsh_coords`, `lhv_polyhedron`, `permissivity_score`

## 📜 License

By contributing code, you agree your work will be released under GPL-3.0.  
Figures and documentation fall under CC BY-NC-ND 4.0.
