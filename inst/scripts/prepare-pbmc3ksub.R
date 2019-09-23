## Process pbmc3k data from TENxPBMCData and subset to highly variable genes
suppressPackageStartupMessages({
    library(TENxPBMCData)
    library(scater)
    library(scran)
    library(BiocSingular)
})

sce <- TENxPBMCData(dataset = "pbmc3k")
colnames(sce) <- paste0("Cell", seq_len(ncol(sce)))
rownames(sce) <- scater::uniquifyFeatureNames(
    ID = rowData(sce)$ENSEMBL_ID,
    names = rowData(sce)$Symbol_TENx
)
counts(sce) <- as.matrix(counts(sce))
sce <- scran::computeSumFactors(sce, min.mean = 0.1)
sce <- scater::normalize(sce)
logcounts(sce) <- as.matrix(logcounts(sce))
new.trend <- scran::makeTechTrend(x = sce)
fit <- scran::trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
fit$trend <- new.trend
dec <- scran::decomposeVar(fit = fit)
top.dec <- dec[order(dec$bio, decreasing = TRUE), ]
keep.genes <- rownames(top.dec[top.dec$bio > 0.01, ])

set.seed(1000)
sce <- scran::denoisePCA(sce, technical = new.trend,
                         subset.row = keep.genes, BSPARAM = IrlbaParam())
reducedDims(sce)$PCA_k2 <- reducedDims(sce)$PCA[, 1:2]

set.seed(1000)
tmp <- scater::runTSNE(sce, use_dimred = "PCA", perplexity = 30)
reducedDims(sce)$TSNE_perp30 <- reducedDims(tmp)$TSNE
set.seed(1000)
tmp <- scater::runTSNE(sce, use_dimred = "PCA", perplexity = 5)
reducedDims(sce)$TSNE_perp5 <- reducedDims(tmp)$TSNE
set.seed(1000)
tmp <- scater::runTSNE(sce, use_dimred = "PCA", perplexity = 15)
reducedDims(sce)$TSNE_perp15 <- reducedDims(tmp)$TSNE
set.seed(1000)
tmp <- scater::runTSNE(sce, use_dimred = "PCA", perplexity = 100)
reducedDims(sce)$TSNE_perp100 <- reducedDims(tmp)$TSNE

set.seed(1000)
tmp <- scater::runUMAP(sce, use_dimred = "PCA", n_neighbors = 15)
reducedDims(sce)$UMAP_nn15 <- reducedDims(tmp)$UMAP
set.seed(1000)
tmp <- scater::runUMAP(sce, use_dimred = "PCA", n_neighbors = 30)
reducedDims(sce)$UMAP_nn30 <- reducedDims(tmp)$UMAP

pbmc3ksub <- sce[keep.genes, ]

## usethis::use_data(pbmc3ksub)
