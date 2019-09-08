#' Evaluate dimension reduction methods
#'
#' Calculate a collection of metrics comparing one or more reduced dimension
#' representations to an underlying high-dimensional data set. The function
#' takes a \code{SingleCellExperiment} object as input, and will use one of the
#' included assays and some or all of the included reduced dimension
#' representations for the evaluation. The "original" distances are calculated
#' from the indicated assay, using the specified variables (or all, if no
#' variables are specified) and a subset of the samples (or all, if nSamples is
#' not specified). These distances are then compared to distances calculated
#' from the specified reduced dimension representations, and several scores are
#' returned. The execution time of the function depends strongly on both the
#' number of retained variables (which affects the distance calculation in the
#' original space) and the number of retained samples. Since subsampling of the
#' columns (via \code{nSamples}) is random, setting the random seed is
#' recommended to obtain reproducible results.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param dimReds A character vector with the names of the reduced dimension
#'   representations from \code{sce} to include in the evaluation.
#' @param assay A character scalar giving the name of the assay from the
#'   \code{sce} to use as the basis for the distance calculations in the
#'   original space.
#' @param features A character vector giving the IDs of the features to use for
#'   distance calculations from the chosen assay. Will be matched to the row
#'   names of the \code{sce}.
#' @param nSamples A numeric scalar, giving the number of columns to subsample
#'   (randomly) from the \code{sce}.
#' @param distNorm A logical scalar, indicating whether the distance vectors in
#'   the original and low-dimensional spaces should be L2 normalized before they
#'   are compared. If set to FALSE, they are instead divided by the square root
#'   of their length, to avoid metrics scaling with the number of retained
#'   samples.
#' @param highDimDistMethod A character scalar defining the distance measure to
#'   use in the original (high-dimensional) space Must be one of "euclidean",
#'   "manhattan", "maximum", "canberra" or "cosine". The distance in the
#'   low-dimensional representation will always be Euclidean.
#' @param kTM An integer vector giving the number of neighbors to use for
#'   trustworthiness and continuity calculations.
#' @param labelColumn A character scalar defining a column of
#'   \code{colData(sce)} to use as the group assignments in the silhouette width
#'   calculations. If not provided, the silhouette width will not be calculated.
#' @param verbose A logical scalar, indicating whether to print out progress
#'   messages.
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @return A \code{data.frame} with values of all evaluation metrics, across
#'   the dimension reduction methods. The following metrics are included:
#'   \itemize{
#'   \item SpearmanCorrDist - The Spearman correlation between the original,
#'   high-dimensional distances and the Euclidean distances in the
#'   low-dimensional representation. Higher values are better.
#'   \item PearsonCorrDist - The Pearson correlation between the original,
#'   high-dimensional distances and the Euclidean distances in the
#'   low-dimensional representation. Higher values are better.
#'   \item KSstatDist - The Kolmogorov-Smirnov statistic comparing the
#'   distribution of distances in the original space and in the low-dimensional
#'   representation. Lower values are better.
#'   \item EuclDistBetweenDists - The Euclidean distance between the vector of
#'   distances in the original, high-dimensional space and those in the
#'   low-dimensional representation. If \code{distNorm} is TRUE, the Euclidean
#'   distance vectors are normalized to have an L2 norm equal to 1 before the
#'   distance between them is calculated. If \code{distNorm} is FALSE, the
#'   distance vectors are instead divided by the square root of their length, to
#'   minimize the influence of the number of samples. Lower values are better.
#'   \item SammonStress - The Sammon stress. If \code{distNorm} is TRUE, the
#'   Euclidean distance vectors are normalized to have an L2 norm equal to 1
#'   before the stress is calculated. If \code{distNorm} is FALSE, the distance
#'   vectors are instead divided by the square root of their length, to minimize
#'   the influence of the number of samples. Lower values are better.
#'   \item Trustworthiness_kNN - The trustworthiness score (Kaski et al. 2003),
#'   using NN nearest neighbors. The trustworthiness indicates to which degree
#'   we can trust that the points placed closest to a given sample in the
#'   low-dimensional representation are really close to the sample also in the
#'   original space. Higher values are better.
#'   \item Continuity_kNN - The continuity score (Kaski et al. 2003), using NN
#'   nearest neighbors. The continuity indicates to which degree we can trust
#'   that the points closest to a given sample in the original space are placed
#'   close to the sample also in the low-dimensional representation. Higher
#'   values are better.
#'   \item MeanJaccard_kNN - The mean Jaccard index (over all samples),
#'   comparing the set of NN nearest neighbors in the original space and those
#'   in the low-dimensional representation. Higher values are better.
#'   \item MeanSilhouette_X - If a \code{labelColumn} X is supplied, the mean
#'   silhouette index across all samples, with the grouping given by this column
#'   and the distances obtained from the low-dimensional representation. Higher
#'   values are better.
#'   \item coRankingQlocal - Q_local as calculated by the coRanking package
#'   (Kraemer and Reichstein 2018, Lee et al. 2009, Chen and Buja 2006). Higher
#'   values are better.
#'   \item coRankingQglobal - Q_global as calculated by the coRanking package
#'   (Kraemer and Reichstein 2018, Lee et al. 2009, Chen and Buja 2006). Higher
#'   values are better.
#'   }
#'
#' @references
#'   Kaski, S., Nikkilä, J., Oja, M., Venna, J., Törönen, P., and
#'   Castrén, E. (2003). Trust- worthiness and metrics in visualizing similarity
#'   of gene expression. BMC Bioinformatics 4:48.
#'
#'   Lee, J.A., Lee, J.A., Verleysen, M. (2009). Quality assessment of
#'   dimensionality reduction: Rankbased criteria. Neurocomputing 72.
#'
#'   Chen, L., Buja, A. (2006). Local Multidimensional Scaling for Nonlinear
#'   Dimension Reduction, Graph Layout and Proximity Analysis.
#'
#'   Kraemer G, Reichstein M, D. MM (2018). dimRed and coRanking-Unifying
#'   Dimensionality Reduction in R. The R Journal 10(1):342-358.
#'
#' @importFrom SingleCellExperiment reducedDims reducedDimNames
#' @importFrom SummarizedExperiment assays colData assayNames
#' @importFrom cluster silhouette
#' @importFrom stats cor ks.test
#' @importFrom dplyr bind_rows
#' @importFrom wordspace dist.matrix
#' @importFrom coRanking coranking LCMC
#' @importFrom methods is
#'
dreval <- function(
    sce, dimReds = NULL, assay = "logcounts",
    features = NULL, nSamples = NULL, distNorm = FALSE,
    highDimDistMethod = "euclidean", kTM = c(10, 100),
    labelColumn = NULL, verbose = FALSE) {

    ## --------------------------------------------------------------------- ##
    ## Check input arguments
    ## --------------------------------------------------------------------- ##
    if (!methods::is(sce, "SingleCellExperiment")) {
        stop("sce must be a SingleCellExperiment object")
    }

    ## If dimReds is NULL, use all dimReds
    if (is.null(dimReds)) {
        dimReds <- SingleCellExperiment::reducedDimNames(sce)
    } else {
        dimReds <- intersect(
            dimReds, SingleCellExperiment::reducedDimNames(sce)
        )
    }
    if (length(dimReds) == 0) (
        stop("No valid reduced dimension representations found")
    )

    if (!(assay %in% SummarizedExperiment::assayNames(sce))) {
        stop("There is no assay named ", assay, " in sce")
    }

    if (!is.logical(distNorm)) {
        stop("distNorm must be a logical value")
    }

    if (!(highDimDistMethod %in% c("euclidean", "manhattan", "maximum",
                                   "canberra", "cosine"))) {
        stop("highDimDistMethod must be one of ",
             "'euclidean', 'manhattan', 'maximum', 'canberra' or 'cosine'")
    }

    if (!is.numeric(kTM)) {
        stop("kTM must be a numeric vector")
    }

    if (!is.logical(verbose)) {
        stop("verbose must be a logical value")
    }

    ## If features is not NULL, subset to provided features
    if (verbose) message("Subsampling rows...")
    if (!is.null(features)) {
        features <- intersect(features, rownames(sce))
        if (length(features) == 0) {
            stop("No features remaining after subsetting. ",
                 "Make sure that the provided features correspond ",
                 "to the row names of sce.")
        }
        sce <- sce[features, ]
    }

    ## If nSamples is not NULL, subsample columns for increased speed
    if (verbose) message("Subsampling columns...")
    if (!is.null(nSamples)) {
        nSamples <- min(nSamples, ncol(sce))
        if (nSamples == 0) {
            stop("nSamples must be >0")
        }
        keep <- sample(seq_len(ncol(sce)), nSamples, replace = FALSE)
        sce <- sce[, keep]
    }

    if (!is.null(labelColumn)) {
        if (!(labelColumn %in%
              colnames(SummarizedExperiment::colData(sce)))) {
            stop("labelColumn must be one of the columns in colData(sce)")
        }
        if (length(unique(SummarizedExperiment::colData(sce)[[labelColumn]])) < 2) {
            stop("labelColumn must correspond to a column ",
                 "with at least two unique values")
        }
        if (length(unique(SummarizedExperiment::colData(sce)[[labelColumn]])) == ncol(sce)) {
            stop("labelColumn must correspond to a column ",
                 "where not all values are different")
        }
    }

    ## --------------------------------------------------------------------- ##
    ## Initialize list to hold results
    ## --------------------------------------------------------------------- ##
    results <- lapply(dimReds, function(m) list())
    names(results) <- dimReds

    ## --------------------------------------------------------------------- ##
    ## Calculate distances and ranks for the high-dimensional data
    ## --------------------------------------------------------------------- ##
    if (verbose) message("Calculating distances in high-dimensional space...")
    mat <- as.matrix(SummarizedExperiment::assays(sce)[[assay]])
    distOriginal <- wordspace::dist.matrix(
        mat, method = highDimDistMethod,
        as.dist = TRUE, byrow = FALSE, convert = FALSE
    )

    ## Get the rank for each sample compared to each other sample. In each
    ## column, the row with the value 1 corresponds to the nearest neighbor,
    ## the one with the value 2 to the second nearest neighbor, etc.
    if (verbose) message("Getting ranks in high-dimensional space...")
    rankOriginal <- apply(
        as.matrix(distOriginal), 2, function(w) order(order(w))
    )

    ## Define normalization constant for distances
    if (distNorm) {
        distNormOriginal <- sqrt(sum(distOriginal^2))
    } else {
        distNormOriginal <- sqrt(length(distOriginal))
    }

    ## --------------------------------------------------------------------- ##
    ## For each dimRed, calculate scores
    ## --------------------------------------------------------------------- ##
    for (dr in dimReds) {
        if (verbose) message(dr, ":")
        dimRedMat <- SingleCellExperiment::reducedDims(sce)[[dr]]

        ## Get dimensionality
        results[[dr]]$dimensionality <- ncol(dimRedMat)

        ## ----------------------------------------------------------------- ##
        ## Calculate distances and ranks for the low-dimensional data
        ## ----------------------------------------------------------------- ##
        if (verbose) message("  Calculating distances in low-dimensional space...")
        distLowDim <- wordspace::dist.matrix(
            as.matrix(dimRedMat), method = "euclidean",
            as.dist = TRUE, byrow = TRUE, convert = FALSE
        )

        ## Get the rank for each sample compared to each other sample. In each
        ## column, the row with the value 1 corresponds to the nearest neighbor,
        ## the one with the value 2 to the second nearest neighbor, etc.
        if (verbose) message("  Getting ranks in low-dimensional space...")
        rankLowDim <- apply(
            as.matrix(distLowDim), 2, function(w) order(order(w))
        )

        if (distNorm) {
            distNormLowDim <- sqrt(sum(distLowDim^2))
        } else {
            distNormLowDim <- sqrt(length(distLowDim))
        }

        ## ----------------------------------------------------------------- ##
        ## Correlation between distance vectors
        ## ----------------------------------------------------------------- ##
        if (verbose) message("  Calculating Spearman correlations...")
        results[[dr]]$SpearmanCorrDist <- stats::cor(
            distOriginal, distLowDim,
            method = "spearman"
        )
        if (verbose) message("  Calculating Pearson correlations...")
        results[[dr]]$PearsonCorrDist <- stats::cor(
            distOriginal, distLowDim,
            method = "pearson"
        )

        ## ----------------------------------------------------------------- ##
        ## Kolmogorov-Smirnov statistic comparing distance distributions
        ## ----------------------------------------------------------------- ##
        if (verbose) message("  Calculating KS statistic")
        results[[dr]]$KSStatDist <- stats::ks.test(
            distOriginal/distNormOriginal,
            distLowDim/distNormLowDim
        )$statistic

        ## ----------------------------------------------------------------- ##
        ## Euclidean distance to original-space (scaled) distances
        ## ----------------------------------------------------------------- ##
        if (verbose) message("  Calculating distances between distances...")
        results[[dr]]$EuclDistBetweenDists <-
            sqrt(sum(((distOriginal/distNormOriginal) -
                          (distLowDim/distNormLowDim))^2))

        ## ----------------------------------------------------------------- ##
        ## Sammon stress
        ## ----------------------------------------------------------------- ##
        if (verbose) message("  Calculating Sammon stress...")
        results[[dr]]$SammonStress <-
            1/sum(distOriginal/distNormOriginal) *
            sum(((distOriginal/distNormOriginal -
                      distLowDim/distNormLowDim) ^ 2)/
                    (distOriginal/distNormOriginal))

        ## ----------------------------------------------------------------- ##
        ## Trustworthiness and continuity
        ## ----------------------------------------------------------------- ##
        for (k in kTM) {
            if (verbose) message("  Calculating trustworthiness, k=", k)
            tm <- calcTrustworthinessFromRank(
                rankOriginal = rankOriginal,
                rankLowDim = rankLowDim,
                kTM = k
            )
            results[[dr]][[paste0("Trustworthiness_k", k)]] <- tm

            if (verbose) message("  Calculating continuity, k=", k)
            ct <- calcContinuityFromRank(
                rankOriginal = rankOriginal,
                rankLowDim = rankLowDim,
                kTM = k
            )
            results[[dr]][[paste0("Continuity_k", k)]] <- ct
        }

        ## ----------------------------------------------------------------- ##
        ## Jaccard index of nearest neighbors
        ## ----------------------------------------------------------------- ##
        for (k in kTM) {
            if (verbose) message("  Calculating Jaccard index ",
                                 "of nearest neighbors, k=", k)
            intrs <- colSums((rankOriginal <= k) * (rankLowDim <= k))
            unin <- colSums(sign((rankOriginal <= k) + (rankLowDim <= k)))
            jaccs <- intrs/unin
            results[[dr]][[paste0("MeanJaccard_k", k)]] <- mean(jaccs)
        }

        ## ----------------------------------------------------------------- ##
        ## Silhouette index relative to known labels
        ## ----------------------------------------------------------------- ##
        if (!is.null(labelColumn)) {
            silh <- cluster::silhouette(
                x = as.numeric(as.factor(
                    SummarizedExperiment::colData(sce)[[labelColumn]]
                )),
                dist = distLowDim
            )
            results[[dr]][[paste0("MeanSilhouette_", labelColumn)]] <-
                summary(silh)$avg.width
        }

        ## ----------------------------------------------------------------- ##
        ## Coranking (with the coRanking package)
        ## ----------------------------------------------------------------- ##
        qt <- coRanking::coranking(Xi = as.matrix(distOriginal),
                                   X = as.matrix(distLowDim),
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
