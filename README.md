# Genetic-Bell (Kaba et al., 2023) 
<p align="left">
  <a href="https://creativecommons.org/licenses/by-nc-nd/4.0/">
    <img src="https://img.shields.io/badge/License-Open_Access-green" alt="">
  </a>
  <a href="https://www.mdpi.com/openaccess">
    <img src="https://img.shields.io/badge/Doi-10.3390/math11244916-blue" alt="">
  </a>
</p>

In complex diseases, the interactions among genes are commonly elucidated through the lens of graphs. Amongst these genes, certain ones form bi-functional modules within the graph, contingent upon their (anti)correlation with a specific functional state, such as susceptibility to a genetic disorder of non-Mendelian traits. Consequently, a disease can be delineated by a finite number of these discernible modules. Within each module, there exist allelic variants that pose a genetic risk, thus qualifying as genetic risk factors. These factors precipitate a permissive state, which if all other modules also align in the same permissive state, can ultimately lead to the onset of the disease in an individual. To gain a deeper insight into the incidence of a disease, it becomes imperative to acquire a comprehensive understanding of the genetic transmission of these factors. In this work, we present a non-linear model for this transmission, drawing inspiration from the classic theory of the Bell experiment. This model aids in elucidating the variances observed in SNP interactions concerning the risk of disease.

```ruby
# Non-linear correlation model of two genes
```
(A) An example of a synthetic biological network modeling a complex disease and clustered by communities, wherein diverse hues signify distinct community affiliations. (B) The correlation model of any two nodes in the same community, affected by potential environmental, stochastic, or genetic variations. A sequence, such as a transcription factor, in cis with gene A, produces a sharp phenotype that interacts with gene B. A variant ω in the gene A promoter (PR) can impede the full phenotypical interaction, leading to a potential disease-permissive state.

<img src="https://github.com/MorillaLab/Genetic-Bell/blob/main/Figure1.png" alt="Non-linear model" width="700">

```ruby
# Genotype versus phenotype co-expression analysis by individual modules and groups
```
(A) Matrix of individual’s GWAS features. (B) Standard scaling methods of a patient matrix, wherein blue and red bars stand for upper and lower threshold of choice respectively. (C) Calculation of the covariance matrix to be subsequently used in data reduction methods. (D) Multidimensional scaling to infer an all-in patient-features matrix. (E) Application of the classic k-means method to separate individuals into two clusters per GO function. Numbers indicate belonging to each cluster

<img src="https://github.com/MorillaLab/Genetic-Bell/blob/main/Figure2.png" alt="Genotype_vs_Phenotype" width="700">

```ruby
# Geometric Bell’s Inequalities of our results
```
The red facets correspond to the 42 dimensional probabilistic CHSH-space of two gene interactions. The green square denotes the potential interaction subspace of (Bax, p21). The LHV region is described by the blue polyhedron. The proven global correlation of (Bax, p21) in functional modules, as applied in modeling complex diseases, places our model in the top-left region between the LHV region and the interactions subspace. Lastly, the non-local and no-signaling regions are indicated at the top-right and bottom-left vertices, respectively.

<img src="https://github.com/MorillaLab/Genetic-Bell/blob/main/Figure7.png" alt="CHSH_space" width="500">


