---
name: Feature / extension request
about: Suggest an extension to the Genetic-Bell model
title: '[FEATURE] '
labels: enhancement
assignees: ''
---

## Motivation
Why would this extend the Genetic-Bell model scientifically?

## Proposed addition
Describe what you'd like: new gene pairs, higher-order interactions, new disease, alternative inequality formulation, improved visualisation, etc.

## Mathematical context (if applicable)
- Proposed formulation or inequality:
- Reference: [Bell 1964 / CHSH 1969 / your extension]

## Biological context
What disease, pathway, or gene network would this apply to?

## Implementation sketch (optional)
```r
# How it would extend the current R API
chsh_3gene <- compute_chsh_3way(data, genes = c("Bax", "p21", "p53"))
```

## References
Relevant papers, databases (GWAS Catalog, GTEx, Hi-C data), or existing tools.
