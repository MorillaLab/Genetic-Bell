---
name: Bug report
about: Report a bug in Genetic-Bell
title: '[BUG] '
labels: bug
assignees: ''
---

## Describe the bug
A clear description of what the bug is.

## Which component fails?
- [ ] Core model (`model/`)
- [ ] EM peak detection (`EM_peaks/`)
- [ ] CHSH space computation
- [ ] Community detection / module analysis
- [ ] Visualisation
- [ ] Other: ___

## Minimal reproducible example
```r
source("model/genetic_bell_model.R")
# Error occurs here:
chsh_coords <- compute_chsh(data, gene_a = "Bax", gene_b = "p21")
```

## Error message
```
Paste full error message here
```

## Data context
- Input data type: [GWAS matrix / synthetic / other]
- Number of genes / SNPs:
- Number of individuals:

## Environment
- OS: [e.g. Ubuntu 22.04 / macOS 14]
- R version: [e.g. 4.3.1] (`R.version`)
- igraph version:
- Key package versions (from `sessionInfo()`)

## Additional context
Any other context about the problem.
