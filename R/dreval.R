#' Evaluate dimension reduction methods
#'
#' Calculate several goodness-of-fit metrics for one or multiple reduced
#' dimension representations. The function takes a \code{SingleCellExperiment}
#' object as input, and will use one of the included assays and some or all of
#' the included reduced dimension representations for the evaluation. The
#' "original" distances are calculated from the indicated assay, using the
#' specified variables (or all, if no variables are specified) and a subset of
#' the samples (or all, if nSamples is not specified). These distances are then
#' compared to distances calculated from the specified reduced dimension
#' representations, and several scores are returned.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param dimReds A character vector giving the names of the reduced dimension
#'   representations from \code{sce} to compare.
#' @param assay A character scalar giving the name of the assay from the
#'   \code{sce} to use as the basis for the comparisons.
#' @param features A character vector giving the IDs of the features to use for
#'   distance calculations from the chosen assay.
#' @param nSamples A numeric scalar, the number of columns to subsample
#'   (randomly) from the \code{sce}.
#' @param distNorm A logical scalar, whether the distance vectors should be L2
#'   normalized before the Sammon stress and the Euclidean distance between them
#'   is calculated. If set to FALSE, they are instead divided by the square root
#'   of their length.
#' @param highDimDistMethod A character scalar defining the distance measure to
#'   use in the original (high-dimensional) data. Must be one of "euclidean",
#'   "manhattan", "maximum", "canberra", "cosine". The distance in the
#'   low-dimensional representation will always be Euclidean.
#' @param kTM An integer vector giving the number of neighbors to use for
#'   trustworthiness and continuity calculations.
#' @param labelColumn A character scalar defining a column of
#'   \code{colData(sce)} to use as the basis for silhouette width calculations.
#' @param verbose A logical scalar, whether to print out progress messages.
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @return A \code{data.frame} with values for all evaluation metrics, across
#'   the dimension reduction methods. The following metrics are included:
#'   \itemize{
#'   \item SpearmanCorrEuclDist - The Spearman correlation between the original,
#'   high-dimensional Euclidean distances and the Euclidean distances in the
#'   low-dimensional representation.
#'   \item PearsonCorrEuclDist - The Pearson correlation between the original,
#'   high-dimensional Euclidean distances and the Euclidean distances in the
#'   low-dimensional representation.
#'   \item EuclDistToEuclDist - The Euclidean distance between the vector of
#'   Euclidean distances in the original, high-dimensional data and those in the
#'   low-dimensional representation. If \code{distNorm} is TRUE, the Euclidean
#'   distance vectors are normalized to have an L2 norm equal to 1 before the
#'   distance between them is calculated. If \code{distNorm} is FALSE, the
#'   distance vectors are instead divided by the square root of their length, to
#'   make the value independent of the number of samples.
#'   \item SammonStress - The Sammon stress. If \code{distNorm} is TRUE, the
#'   Euclidean distance vectors are normalized to have an L2 norm equal to 1
#'   before the stress is calculated. If \code{distNorm} is FALSE, the distance
#'   vectors are instead divided by the square root of their length, to make the
#'   value independent of the number of samples.
#'   \item TrustworthinessEuclDist_kNN - The trustworthiness score (Kaski et al.
#'   2003), using NN nearest neighbors. The trustworthiness indicates to which
#'   degree we can trust that the points placed closest to a given sample in the
#'   low-dimensional representation are really close to the sample also in the
#'   original data set.
#'   \item ContinuityEuclDist_kNN - The continuity score (Kaski et al. 2003),
#'   using NN nearest neighbors. The continuity indicates to which degree we can
#'   trust that the points closest to a given sample in the original data set
#'   are placed close to the sample also in the low-dimensional representation.
#'   \item MeanJaccard_kNN - The mean Jaccard index (over all samples),
#'   comparing the set of NN nearest neighbors in the original data and those in
#'   the low-dimensional representation.
#'   \item MeanSilhouette_X - If a \code{labelColumn} is supplied, the mean
#'   silhouette index across all samples, with the grouping given by this column
#'   and the distances obtained from the low-dimensional representation.
#'   \item coRankingQlocal - Q_local as calculated by the coRanking package.
#'   \item coRankingQglobal - Q_global as calculated by the coRanking package.
#'   }
#'
#' @references Kaski, S., Nikkilä, J., Oja, M., Venna, J., Törönen, P., and
#'   Castrén, E. (2003). Trust- worthiness and metrics in visualizing similarity
#'   of gene expression. BMC Bioinformatics 4:48.
#'
#' @importFrom SingleCellExperiment reducedDims reducedDimNames
#' @importFrom SummarizedExperiment assays colData
#' @importFrom cluster silhouette
#' @importFrom stats cor
#' @importFrom dplyr bind_rows
#' @importFrom wordspace dist.matrix
#' @importFrom coRanking coranking LCMC
#'
dreval <- function(
    sce, dimReds = NULL, assay = "logcounts",
    features = NULL, nSamples = NULL, distNorm = FALSE,
    highDimDistMethod = "euclidean", kTM = c(10, 100),
    labelColumn = NULL, verbose = FALSE) {

    ## Initialize list to hold results
    results <- lapply(dimReds, function(m) list())
    names(results) <- dimReds

    ## If nSamples is not NULL, subsample columns for increased speed
    if (verbose) message("Subsampling columns...")
    if (!is.null(nSamples)) {
        nSamples <- min(nSamples, ncol(sce))
        keep <- sample(seq_len(ncol(sce)), nSamples, replace = FALSE)
        sce <- sce[, keep]
    }

    ## If features is not NULL, subset to provided features
    if (verbose) message("Subsampling rows...")
    if (!is.null(features)) {
        features <- intersect(features, rownames(sce))
        sce <- sce[features, ]
    }

    ## If dimReds is NULL, use all dimReds
    if (is.null(dimReds)) {
        dimReds <- SingleCellExperiment::reducedDimNames(sce)
    } else {
        dimReds <- intersect(
            dimReds, SingleCellExperiment::reducedDimNames(sce)
        )
    }

    if (is.null(assay)) {
        assay <- 1
    }

    ## Calculate Euclidean distances for the original data
    if (verbose) message("Calculating original Euclidean distances...")
    mat <- as.matrix(SummarizedExperiment::assays(sce)[[assay]])
    euclDistOriginal <- wordspace::dist.matrix(mat, method = highDimDistMethod,
                                               as.dist = TRUE, byrow = FALSE,
                                               convert = FALSE)
    # euclDistOriginal <- stats::dist(t(mat))

    ## Get the rank for each sample compared to each other sample In each
    ## column, the row with the value 1 corresponds to the nearest neighbor,
    ## the one with the value 2 to the second nearest neighbor, etc.
    if (verbose) message("Getting ranks for original data...")
    euclRankOriginal <- apply(
        as.matrix(euclDistOriginal), 2, function(w) order(order(w))
    )

    ## For each dimRed, calculate scores
    for (dr in dimReds) {
        if (verbose) message(dr, ":")
        dimRedMat <- SingleCellExperiment::reducedDims(sce)[[dr]]

        ## Get dimension
        results[[dr]]$dimensionality <- ncol(dimRedMat)

        ## Euclidean distances
        if (verbose) message("  Calculating Euclidean distances...")
        euclDistLowDim <- wordspace::dist.matrix(as.matrix(dimRedMat), method = "euclidean",
                                                 as.dist = TRUE, byrow = TRUE,
                                                 convert = FALSE)
        # euclDistLowDim <- stats::dist(dimRedMat)

        ## Correlation with original-space distances
        if (verbose) message("  Calculating Spearman correlations...")
        results[[dr]]$SpearmanCorrEuclDist <- stats::cor(
            euclDistOriginal, euclDistLowDim,
            method = "spearman"
        )
        if (verbose) message("  Calculating Pearson correlations...")
        results[[dr]]$PearsonCorrEuclDist <- stats::cor(
            euclDistOriginal, euclDistLowDim,
            method = "pearson"
        )

        ## Euclidean distance to original-space (scaled) distances
        if (verbose) message("  Calculating distances between distances...")
        if (distNorm) {
            distNormOriginal <- sqrt(sum(euclDistOriginal^2))
            distNormLowDim <- sqrt(sum(euclDistLowDim^2))
        } else {
            distNormOriginal <- sqrt(length(euclDistOriginal))
            distNormLowDim <- sqrt(length(euclDistLowDim))
        }
        results[[dr]]$EuclDistToEuclDist <-
            sqrt(sum(((euclDistOriginal/distNormOriginal) -
                          (euclDistLowDim/distNormLowDim))^2))

        ## Sammon stress
        if (verbose) message("  Calculating Sammon stress...")
        results[[dr]]$SammonStress <-
            1/sum(euclDistOriginal/distNormOriginal) *
            sum(((euclDistOriginal/distNormOriginal -
                      euclDistLowDim/distNormLowDim) ^ 2)/
                    (euclDistOriginal/distNormOriginal))

        ## Trustworthiness and continuity
        if (verbose) message("  Calculating ranks from low-dim data...")
        euclRankLowDim <- apply(as.matrix(euclDistLowDim), 2,
                                function(w) order(order(w)))
        for (k in kTM) {
            if (verbose) message("  Calculating trustworthiness, k=", k)
            tm <- calcTrustworthinessFromRank(
                rankOriginal = euclRankOriginal,
                rankLowDim = euclRankLowDim,
                kTM = k
            )
            results[[dr]][[paste0("TrustworthinessEuclDist_k", k)]] <- tm

            if (verbose) message("  Calculating continuity, k=", k)
            ct <- calcContinuityFromRank(
                rankOriginal = euclRankOriginal,
                rankLowDim = euclRankLowDim,
                kTM = k
            )
            results[[dr]][[paste0("ContinuityEuclDist_k", k)]] <- ct
        }

        ## Jaccard index of nearest neighbors
        for (k in kTM) {
            if (verbose) message("  Calculating Jaccard index of nearest neighbors, k=", k)
            intrs <- colSums((euclRankOriginal <= k) * (euclRankLowDim <= k))
            unin <- colSums(sign((euclRankOriginal <= k) + (euclRankLowDim <= k)))
            jaccs <- intrs/unin
            results[[dr]][[paste0("MeanJaccard_k", k)]] <- mean(jaccs)
        }

        ## Silhouette index relative to known labels
        if (!is.null(labelColumn) && labelColumn %in%
            colnames(SummarizedExperiment::colData(sce)) &&
            length(unique(SummarizedExperiment::colData(sce)[[labelColumn]])) >= 2) {
            silh <- cluster::silhouette(
                x = as.numeric(as.factor(SummarizedExperiment::colData(sce)[[labelColumn]])),
                dist = euclDistLowDim
            )
            results[[dr]][[paste0("MeanSilhouette_", labelColumn)]] <- summary(silh)$avg.width
        }

        ## Coranking (with the coRanking package)
        qt <- coRanking::coranking(Xi = as.matrix(euclDistOriginal),
                                   X = as.matrix(euclDistLowDim),
                                   input = "dist")
        lcmc <- coRanking::LCMC(qt)
        Kmax <- which.max(lcmc)
        qlocal <- mean(lcmc[seq(from = 1, to = Kmax, by = 1)])
        qglobal <- mean(lcmc[seq(from = Kmax + 1, to = length(lcmc), by = 1)])
        results[[dr]][[paste0("coRankingQlocal")]] <- qlocal
        results[[dr]][[paste0("coRankingQglobal")]] <- qglobal
    }

    results <- do.call(dplyr::bind_rows, lapply(names(results), function(nm) {
        data.frame(
            Method = nm, as.data.frame(results[[nm]]),
            stringsAsFactors = FALSE
        )
    }))
    results
}
