context("plotting")

test_that("plotting works", {
    data(pbmc3ksub)
    set.seed(1)
    sce <- pbmc3ksub[1:100, 1:100]
    dres <- dreval(
        sce = sce, dimReds = c("PCA", "PCA_k2"),
        refType = "assay", refAssay = "logcounts",
        features = NULL, nSamples = NULL, distNorm = "l2",
        refDistMethod = "euclidean", kTM = 5,
        labelColumn = NULL, verbose = FALSE
    )

    gg <- plotRankSummary(dres$scores, sortBars = "decreasing")
    expect_equal(gg$data$score[gg$data$Method == "PCA"],
                 rep(2, 10))
    expect_equal(gg$data$score[gg$data$Method == "PCA_k2"],
                 rep(1, 10))
    expect_equal(gg$data$Method, factor(rep(c("PCA", "PCA_k2"), each = 10),
                                        levels = c("PCA", "PCA_k2")))

    gg <- plotRankSummary(dres$scores, sortBars = "increasing")
    expect_equal(gg$data$score[gg$data$Method == "PCA"],
                 rep(2, 10))
    expect_equal(gg$data$score[gg$data$Method == "PCA_k2"],
                 rep(1, 10))
    expect_equal(gg$data$Method, factor(rep(c("PCA_k2", "PCA"), each = 10),
                                        levels = c("PCA_k2", "PCA")))

    gg <- plotRankSummary(dres$scores, sortBars = "decreasing",
                          metrics = "global")
    expect_equal(gg$data$score[gg$data$Method == "PCA"],
                 rep(2, 6))
    expect_equal(gg$data$score[gg$data$Method == "PCA_k2"],
                 rep(1, 6))
    expect_equal(gg$data$Method, factor(rep(c("PCA", "PCA_k2"), each = 6),
                                        levels = c("PCA", "PCA_k2")))
    expect_equal(sort(unique(as.character(gg$data$metric))),
                 sort(c("coRankingQglobal", "EuclDistBetweenDists", "KSStatDist",
                        "PearsonCorrDist", "SammonStress", "SpearmanCorrDist")))

    gg <- plotRankSummary(dres$scores, sortBars = "decreasing",
                          metrics = "local")
    expect_equal(gg$data$score[gg$data$Method == "PCA"],
                 rep(2, 4))
    expect_equal(gg$data$score[gg$data$Method == "PCA_k2"],
                 rep(1, 4))
    expect_equal(gg$data$Method, factor(rep(c("PCA", "PCA_k2"), each = 4),
                                        levels = c("PCA", "PCA_k2")))
    expect_equal(sort(unique(as.character(gg$data$metric))),
                 sort(c("Continuity_k5", "Trustworthiness_k5",
                        "coRankingQlocal", "MeanJaccard_k5")))
})
