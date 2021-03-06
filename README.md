## dreval
[![Codecov.io coverage status](https://codecov.io/github/csoneson/dreval/coverage.svg?branch=master)](https://codecov.io/github/csoneson/dreval)
[![R build status](https://github.com/csoneson/dreval/workflows/R-CMD-check/badge.svg)](https://github.com/csoneson/dreval)

`dreval` is an R package aimed at evaluation and comparison of reduced dimension
representations of high-dimensional data. Given one or more reduced dimension
representations, and a "reference" representation (which can be the original,
high-dimensional representation or a baseline low-dimensional one), `dreval`
will calculate a collection of metrics quantifying how well each of the
evaluated representations recapitulates the structure of the observations in the
reference representation.

### Installation

To install `dreval`, you need the `remotes` (or `devtools`) R package, which can
be installed from CRAN. The following commands installs first `remotes`, then
`dreval`.

```
install.packages("remotes")
remotes::install_github("csoneson/dreval")
```

### Application

The input to `dreval` is a
[`SingleCellExperiment`](https://bioconductor.org/packages/SingleCellExperiment/)
object, containing one or more assays and one or more reduced dimension
representation. By default, the `logcounts` assay will be used as the reference
representation, against which each of the provided reduced dimension
representations will be evaluated. However, any other assay or reduced dimension
representation can be used as the reference data, by setting the arguments to
the `dreval()` function accordingly.

The package contains a small example single-cell RNA-seq data set with
measurements for approximately 1,800 highly variable genes across 2,700 PBMCs.
The object contains eight reduced dimension representations: 25-dimensional PCA,
2-dimensional PCA, and 2-dimensional t-SNE and 2-dimensional UMAP
representations generated with different values of the perplexity/number of
nearest neighbors. We use the `dreval()` function to evaluate how well each of
these retain the structure of the cells based on the `logcounts` assay.

```
data(pbmc3ksub)
dre <- dreval(sce = pbmc3ksub, refType = "assay", 
              refAssay = "logcounts", nSamples = 1000, kTM = 50)
```

For detailed information about the arguments to `dreval()` we refer to the
help page of the function:

```
?dreval
```

The output of `dreval()` is a list with two elements, named `scores` and
`plots`. The `scores` element is a data.frame with all the calculated evaluation
scores for each of the reduced dimension representations, while the `plots`
element is a list of diagnostic plots.

The `plotRankSummary()` function can be used to aggregate the information across
all evaluation metrics. Each reduced dimension representation will be assigned a
rank for each metric, and the sum of these ranks across all metrics, as well as
the contribution from each metric, is visualized by the function. Metrics aimed
at evaluating the preservation of global structure are colored blue, while those
aimed at evaluating the preservation of local structure are colored red.

```
plotRankSummary(dre$scores)
```

### Related material

- A python framework for reduced dimension representation evaluation was proposed by [Heiser and Lau (2019)](https://www.biorxiv.org/content/10.1101/684340v1.abstract). Code is available on [GitHub](https://github.com/KenLauLab/DR-structure-preservation). This study proposed to use the Mantel correlation between distance matrices, the Earth Mover's Distance between distance distributions, and the percentage of total binary elements of the KNN matrix that are conserved as evaluation metrics. 
- The [dimRed](https://cran.r-project.org/web/packages/dimRed/index.html) R package implements a collection of quality scores for reduced dimension representations, including Q\_local, Q\_global (based on the co-ranking matrix) and various correlation measures. 
- The [quadra](https://github.com/jlmelville/quadra) R package implements quality scores for reduced dimension representations based on preservation of neighborhoods and agreement with known labels. 

