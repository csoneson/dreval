## dreval

`dreval` is an R package for evaluation and comparison of reduced dimension representations, in terms of how well they represent the structure of the underlying high-dimensional data set.

### Installation

To install `dreval`, you need the `remotes` R package, which can be installed from CRAN.

```
install.packages("remotes")
remotes::install_github("csoneson/dreval")
```

### Application

The input to `dreval` is a [`SingleCellExperiment`](https://bioconductor.org/packages/SingleCellExperiment/) object, containing one or more assays and one or more reduced dimension representation. 

### Related material

- A python framework for reduced dimension representation evaluation was proposed by [Heiser and Lau (2019)](https://www.biorxiv.org/content/10.1101/684340v1.abstract). Code is available on [GitHub](https://github.com/KenLauLab/DR-structure-preservation). This study proposed to use the Mantel correlation between distance matrices, the Earth Mover's Distance between distance distributions, and the percentage of total binary elements of the KNN matrix that are conserved as evaluation metrics. 
- The [dimRed]() R package implements a collection of quality scores for reduced dimension representations, including Q\_local, Q\_global (based on the co-ranking matrix) and various correlation measures. 