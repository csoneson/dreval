context("score calculations")

test_that("scores are calculated correctly", {
    data(pbmc3ksub)
    set.seed(1)
    sce <- pbmc3ksub[1:100, 1:100]
    dres <- dreval(
        sce = sce, dimReds = "PCA", refType = "assay", refAssay = "logcounts",
        features = NULL, nSamples = NULL, distNorm = "l2",
        refDistMethod = "euclidean", kTM = 5,
        labelColumn = NULL, verbose = FALSE
    )

    expect_equal(nrow(dres$scores), 1)

    ## recalculate scores
    dists_orig <- wordspace::dist.matrix(
        SingleCellExperiment::logcounts(sce), method = "euclidean",
        byrow = FALSE, as.dist = TRUE
    )
    dists_orig <- dists_orig/sqrt(sum(dists_orig ^ 2))
    ranks_orig <- apply(as.matrix(dists_orig), 2, function(w) order(order(w)))
    expect_equal(diag(ranks_orig), rep(1, ncol(ranks_orig)))

    dists_lowdim <- wordspace::dist.matrix(
        SingleCellExperiment::reducedDim(sce, "PCA"), method = "euclidean",
        byrow = TRUE, as.dist = TRUE
    )
    dists_lowdim <- dists_lowdim/sqrt(sum(dists_lowdim ^ 2))
    ranks_lowdim <- apply(as.matrix(dists_lowdim), 2, function(w) order(order(w)))
    expect_equal(diag(ranks_lowdim), rep(1, ncol(ranks_lowdim)))

    expect_equal(dim(ranks_orig), dim(ranks_lowdim))
    expect_equal(length(dists_orig), length(dists_lowdim))

    N <- ncol(ranks_orig)
    k <- 5

    expect_equal(
        dres$scores$SpearmanCorrDist,
        cor(dists_orig, dists_lowdim, method = "spearman")
    )
    expect_equal(
        dres$scores$PearsonCorrDist,
        cor(dists_orig, dists_lowdim, method = "pearson")
    )
    expect_equivalent(
        dres$scores$KSStatDist,
        ks.test(dists_orig, dists_lowdim)$statistic
    )
    expect_equal(
        dres$scores$EuclDistBetweenDists,
        sqrt(sum((dists_orig - dists_lowdim) ^ 2))
    )
    expect_equal(
        dres$scores$SammonStress,
        1/sum(dists_orig) * sum(((dists_orig - dists_lowdim) ^ 2)/dists_orig)
    )
    expect_equal(
        dres$scores$Trustworthiness_k5,
        1 - 2/(N * k * (2 * N - 3 * k - 1)) *
            sum(sapply(1:ncol(ranks_orig), function(i) {
                (ranks_orig[, i] - k) * (ranks_lowdim[, i] <= k) *
                    (ranks_orig[, i] > k)
            }))
    )
    expect_equal(
        dres$scores$Trustworthiness_k5,
        calcTrustworthinessFromDist(dists_orig, dists_lowdim, kTM = k)
    )
    expect_equal(
        dres$scores$Continuity_k5,
        1 - 2/(N * k * (2 * N - 3 * k - 1)) *
            sum(sapply(1:ncol(ranks_lowdim), function(i) {
                (ranks_lowdim[, i] - k) * (ranks_orig[, i] <= k) *
                    (ranks_lowdim[, i] > k)
            }))
    )
    expect_equal(
        dres$scores$Continuity_k5,
        calcContinuityFromDist(dists_orig, dists_lowdim, kTM = k)
    )
    expect_equal(
        dres$scores$MeanJaccard_k5,
        mean(sapply(1:ncol(ranks_orig), function(i) {
            length(intersect(which(ranks_orig[, i] <= k), which(ranks_lowdim[, i] <= k))) /
                length(union(which(ranks_orig[, i] <= k), which(ranks_lowdim[, i] <= k)))
        }))
    )

    qt <- coRanking::coranking(Xi = as.matrix(dists_orig),
                               X = as.matrix(dists_lowdim),
                               input = "dist")
    lcmc <- coRanking::LCMC(qt)
    Kmax <- which.max(lcmc)
    expect_equal(
        dres$scores$coRankingQglobal,
        mean(lcmc[seq(from = Kmax + 1, to = length(lcmc), by = 1)])
    )
    expect_equal(
        dres$scores$coRankingQlocal,
        mean(lcmc[seq(from = 1, to = Kmax, by = 1)])
    )
})
