#' Evaluate structure preservation in reduced dimension representations
#'
#' Calculates a collection of metrics comparing one or more reduced dimension
#' representations to a reference representation. The function takes a
#' \code{SingleCellExperiment} object as input. The reference representation can
#' be either one of the included assays or one of the reduced dimension
#' representations. If an assay is used, reference distances can be calculated
#' based on all or a subset of the features (rows). These distances are then
#' compared to distances calculated from the specified reduced dimension
#' representations, and several scores are returned. The execution time of the
#' function depends strongly on both the number of retained variables (which
#' affects the distance calculation in the reference space) and the number of
#' samples that are randomly selected to use as the basis for the comparison.
#' Since subsampling of the columns (via the \code{nSamples} argument) is
#' random, setting the random seed is recommended to obtain reproducible
#' results.
#'
#' The following metrics are calculated:
#'   \itemize{
#'   \item SpearmanCorrDist - The Spearman correlation between the reference
#'   distances and the Euclidean distances in the low-dimensional
#'   representation. Ranges from -1 to 1, higher values are better.
#'   \item PearsonCorrDist - The Pearson correlation between the reference
#'   distances and the Euclidean distances in the low-dimensional
#'   representation. Ranges from -1 to 1, higher values are better.
#'   \item KSstatDist - The Kolmogorov-Smirnov statistic comparing the
#'   distribution of distances in the reference space and in the low-dimensional
#'   representation. Ranges from 0 to 1, lower values are better.
#'   \item EuclDistBetweenDists - The Euclidean distance between the vector of
#'   distances in the reference space and those in the low-dimensional
#'   representation. Depending on the value of \code{distNorm}, distances are
#'   scaled before they are compared. Lower values are better.
#'   \item SammonStress - The Sammon stress (Sammon 1969). Depending on the
#'   value of \code{distNorm}, distances are scaled before they are compared.
#'   Lower values are better.
#'   \item Trustworthiness_kNN - The trustworthiness score (Venna & Kaski 2001),
#'   using NN nearest neighbors. The trustworthiness indicates to which degree
#'   we can trust that the points placed closest to a given sample in the
#'   low-dimensional representation are really close to the sample also in the
#'   reference space. Ranges from 0 to 1, higher values are better.
#'   \item Continuity_kNN - The continuity score (Venna & Kaski 2001), using NN
#'   nearest neighbors. The continuity indicates to which degree we can trust
#'   that the points closest to a given sample in the reference space are placed
#'   close to the sample also in the low-dimensional representation. Ranges from
#'   0 to 1, higher values are better.
#'   \item MeanJaccard_kNN - The mean Jaccard index (over all samples),
#'   comparing the set of NN nearest neighbors in the reference space and those
#'   in the low-dimensional representation. Ranges from 0 to 1, higher values
#'   are better.
#'   \item MeanSilhouette_X - If a \code{labelColumn} X is supplied, the mean
#'   silhouette score (Rousseeuw 1987) across all samples, with the grouping
#'   given by this column and the distances obtained from the low-dimensional
#'   representation. Ranges from -1 to 1, higher values are better.
#'   \item coRankingQlocal - Q_local, defined as the average LCMC over the
#'   values to the left of the maximum, following the dimRed/coRanking package
#'   implementations (Kraemer et al 2018, Lee and Verleysen 2009, Chen and Buja
#'   2009). Measures the preservation of local distances, higher values are
#'   better.
#'   \item coRankingQglobal - Q_global, defined as the average LCMC over the
#'   values to the right of the maximum, following the dimRed/coRanking package
#'   implementations (Kraemer et al 2018, Lee and Verleysen 2009, Chen and Buja
#'   2009). Measures the preservation of global distances, higher values are
#'   better.
#'   }
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param dimReds A character vector with the names of the reduced dimension
#'   representations from \code{sce} to include in the evaluation. If
#'   \code{NULL}, all reduced dimension representations are included.
#' @param refType A character scalar, either "assay" or "dimred", specifying
#'   whether to use an assay or a reduced dimension representation of \code{sce}
#'   as the reference data source.
#' @param refAssay A character scalar giving the name of the assay from
#'   \code{sce} to use as the basis for the distance calculations in the
#'   reference space, if \code{refType} if \code{"assay"}.
#' @param refDimRed A character scalar specifying the reduced dimension
#'   representation to use as the reference data representation if
#'   \code{refType} is \code{"dimred"}.
#' @param features A character vector giving the IDs of the features to use for
#'   distance calculations from the chosen assay. Will be matched to the row
#'   names of \code{sce}.
#' @param nSamples A numeric scalar, giving the number of columns to subsample
#'   (randomly) from \code{sce}.
#' @param distNorm A character scalar, indicating how the distance vectors in
#'   the reference and low-dimensional spaces should be normalized before they
#'   are compared. If set to "l2", the vectors are L2 normalized, if set to
#'   "median" they are divided by the median value times the square root of
#'   their length, and if set to any other value they are divided by the square
#'   root of their length, to avoid metrics scaling with the number of retained
#'   samples.
#' @param refDistMethod A character scalar defining the distance measure to use
#'   in the reference space. Must be one of "euclidean", "manhattan", "maximum",
#'   "canberra" or "cosine". The distance in the low-dimensional representation
#'   will always be Euclidean.
#' @param kTM An integer vector giving the number of neighbors to use for
#'   trustworthiness, continuity and Jaccard index calculations.
#' @param labelColumn A character scalar defining a column of
#'   \code{colData(sce)} to use as the group assignments in the silhouette width
#'   calculations. If not provided, the silhouette widths are not calculated.
#' @param verbose A logical scalar, indicating whether to print out progress
#'   messages.
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @return A list with two elements:
#' \itemize{
#' \item scores - A \code{data.frame} with values of all evaluation metrics,
#' across the dimension reduction methods. In addition to the metrics, it
#' contains the dimensionality of the respective reduced dimension
#' representations, and the value of K giving the highest value of LCMC (used
#' for the calculations of Qlocal and Qglobal, see Kraemer et al 2018, Lee and
#' Verleysen 2009, Chen and Buja 2009).
#' \item plots - A list of ggplot objects, representing diagnostic plots.
#' }
#'
#' @references
#'   Venna J., Kaski S. (2001). Neighborhood preservation in nonlinear
#'   projection methods: An experimental study. In Dorffner G., Bischof H.,
#'   Hornik K., editors, Proceedings of ICANN 2001, pp 485–491. Springer,
#'   Berlin.
#'
#'   Lee J.A., Verleysen M. (2009). Quality assessment of
#'   dimensionality reduction: Rank-based criteria. Neurocomputing 72
#'   (7-9):1431-1443.
#'
#'   Chen L., Buja A. (2009). Local multidimensional scaling for nonlinear
#'   dimension reduction, graph drawing, and proximity analysis. Journal of the
#'   American Statistical Association 104:209-219.
#'
#'   Kraemer G., Reichstein M., Mahecha M.D. (2018). dimRed and coRanking -
#'   Unifying dimensionality reduction in R. The R Journal 10 (1):342-358.
#'
#'   Sammon J.W. Jr (1969). A nonlinear mapping for data structure analysis.
#'   IEEE Transactions on Computers C18(5):401-409.
#'
#'   Rousseeuw, P.J. (1987). Silhouettes: A graphical aid to the interpretation
#'   and validation of cluster analysis. Journal of Computational and Applied
#'   Mathematics 20:53-65.
#'
#' @importFrom SingleCellExperiment reducedDims reducedDimNames
#' @importFrom SummarizedExperiment assays colData assayNames
#' @importFrom cluster silhouette
#' @importFrom stats cor ks.test median
#' @importFrom dplyr bind_rows
#' @importFrom wordspace dist.matrix
#' @importFrom coRanking coranking LCMC
#' @importFrom methods is
#' @importFrom ggplot2 ggplot aes geom_hex theme_bw labs scale_fill_gradient
#'   geom_abline stat_ecdf
#' @importFrom tidyr gather
#'
#' @examples
#' data(pbmc3ksub)
#' dre <- dreval(sce = pbmc3ksub, nSamples = 150)
#'
dreval <- function(
    sce, dimReds = NULL, refType = "assay",
    refAssay = "logcounts", refDimRed = NULL,
    features = NULL, nSamples = NULL, distNorm = "none",
    refDistMethod = "euclidean", kTM = c(10, 100),
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
    ## If refType is dimred, don't compare the reference dimred to itself
    if (refType == "dimred") {
        dimReds <- setdiff(dimReds, refDimRed)
    }
    if (length(dimReds) == 0) (
        stop("No valid reduced dimension representations found")
    )

    if (!(refType %in% c("assay", "dimred"))) {
        stop("'refType' must be either 'assay' or 'dimred'")
    }

    if (refType == "assay" &&
        !(refAssay %in% SummarizedExperiment::assayNames(sce))) {
        stop("There is no assay named ", refAssay, " in sce")
    }

    if (refType == "dimred" &&
        !(refDimRed %in% SingleCellExperiment::reducedDimNames(sce))) {
        stop("There is no reduced dimension representation named ",
             refDimRed, " in sce")
    }

    if (!is.character(distNorm)) {
        stop("distNorm must be a character string")
    }

    if (!(refDistMethod %in% c("euclidean", "manhattan", "maximum",
                               "canberra", "cosine"))) {
        stop("refDistMethod must be one of ",
             "'euclidean', 'manhattan', 'maximum', 'canberra' or 'cosine'")
    }

    if (!is.numeric(kTM)) {
        stop("kTM must be a numeric vector")
    }

    if (!is.logical(verbose)) {
        stop("verbose must be a logical value")
    }

    ## If features is not NULL, subset to provided features
    if (refType == "assay" && !is.null(features)) {
        if (verbose) message("Subsampling rows...")
        features <- intersect(features, rownames(sce))
        if (length(features) == 0) {
            stop("No features remaining after subsetting. ",
                 "Make sure that the provided features correspond ",
                 "to the row names of sce.")
        }
        sce <- sce[features, ]
    }

    ## If nSamples is not NULL, subsample columns for increased speed
    if (!is.null(nSamples)) {
        if (verbose) message("Subsampling columns...")
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
        tmp <- SummarizedExperiment::colData(sce)[[labelColumn]]
        if (length(unique(tmp)) < 2) {
            stop("labelColumn must correspond to a column ",
                 "with at least two unique values")
        }
        if (length(unique(tmp)) == ncol(sce)) {
            stop("labelColumn must correspond to a column ",
                 "where not all values are different")
        }
    }

    ## --------------------------------------------------------------------- ##
    ## Initialize list to hold results
    ## --------------------------------------------------------------------- ##
    results <- lapply(dimReds, function(m) list())
    names(results) <- dimReds

    plots <- list(disthex = list())

    ## --------------------------------------------------------------------- ##
    ## Calculate distances and ranks for the high-dimensional data
    ## --------------------------------------------------------------------- ##
    if (verbose) message("Calculating distances in high-dimensional space...")
    if (refType == "assay") {
        mat <- as.matrix(SummarizedExperiment::assays(sce)[[refAssay]])
        distReference <- wordspace::dist.matrix(
            mat, method = refDistMethod,
            as.dist = TRUE, byrow = FALSE, convert = FALSE
        )
    } else if (refType == "dimred") {
        mat <- SingleCellExperiment::reducedDims(sce)[[refDimRed]]
        distReference <- wordspace::dist.matrix(
            mat, method = refDistMethod,
            as.dist = TRUE, byrow = TRUE, convert = FALSE
        )
    } else {
        stop("Invalid 'refType'")
    }

    ## Get the rank for each sample compared to each other sample. In each
    ## column, the row with the value 1 corresponds to the nearest neighbor,
    ## the one with the value 2 to the second nearest neighbor, etc.
    if (verbose) message("Getting ranks in high-dimensional space...")
    rankReference <- apply(
        as.matrix(distReference), 2, function(w) order(order(w))
    )

    ## Define normalization constant for distances
    if (distNorm == "l2") {
        distNormReference <- sqrt(sum(distReference^2))
    } else if (distNorm == "median") {
        distNormReference <- stats::median(distReference) *
            sqrt(length(distReference))
    } else {
        distNormReference <- sqrt(length(distReference))
    }

    distdf <- data.frame(reference = c(distReference/distNormReference))

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
        if (verbose)
            message("  Calculating distances in low-dimensional space...")
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

        if (distNorm == "l2") {
            distNormLowDim <- sqrt(sum(distLowDim^2))
        } else if (distNorm == "median") {
            distNormLowDim <- stats::median(distLowDim) *
                sqrt(length(distLowDim))
        } else {
            distNormLowDim <- sqrt(length(distLowDim))
        }

        distdf[[dr]] <- c(distLowDim/distNormLowDim)

        ## ----------------------------------------------------------------- ##
        ## Correlation between distance vectors
        ## ----------------------------------------------------------------- ##
        if (verbose) message("  Calculating Spearman correlations...")
        results[[dr]]$SpearmanCorrDist <- stats::cor(
            distReference, distLowDim,
            method = "spearman"
        )
        if (verbose) message("  Calculating Pearson correlations...")
        results[[dr]]$PearsonCorrDist <- stats::cor(
            distReference, distLowDim,
            method = "pearson"
        )

        plots$disthex[[dr]] <-
            ggplot2::ggplot(
                data.frame(reference = c(distReference/distNormReference),
                           lowdim = c(distLowDim/distNormLowDim)),
                ggplot2::aes(x = reference, y = lowdim)
            ) +
            ggplot2::geom_hex(bins = 100, ggplot2::aes(fill = stat(density))) +
            ggplot2::scale_fill_gradient(
                name = "", low = "bisque2", high = "darkblue"
            ) +
            ggplot2::theme_bw() +
            ggplot2::labs(
                title = dr, x = "Scaled distance in reference space",
                y = "Scaled distance in low-dimensional space"
            ) +
            ggplot2::geom_abline(slope = 1, intercept = 0)

        ## ----------------------------------------------------------------- ##
        ## Kolmogorov-Smirnov statistic comparing distance distributions
        ## ----------------------------------------------------------------- ##
        if (verbose) message("  Calculating KS statistic")
        results[[dr]]$KSStatDist <- stats::ks.test(
            distReference/distNormReference,
            distLowDim/distNormLowDim
        )$statistic

        ## ----------------------------------------------------------------- ##
        ## Euclidean distance to reference-space (scaled) distances
        ## ----------------------------------------------------------------- ##
        if (verbose) message("  Calculating distances between distances...")
        results[[dr]]$EuclDistBetweenDists <-
            sqrt(sum(((distReference/distNormReference) -
                          (distLowDim/distNormLowDim))^2))

        ## ----------------------------------------------------------------- ##
        ## Sammon stress
        ## ----------------------------------------------------------------- ##
        if (verbose) message("  Calculating Sammon stress...")
        results[[dr]]$SammonStress <-
            1/sum(distReference/distNormReference) *
            sum(((distReference/distNormReference -
                      distLowDim/distNormLowDim) ^ 2) /
                    (distReference/distNormReference))

        ## ----------------------------------------------------------------- ##
        ## Trustworthiness and continuity
        ## ----------------------------------------------------------------- ##
        for (k in kTM) {
            if (verbose) message("  Calculating trustworthiness, k=", k)
            tm <- calcTrustworthinessFromRank(
                rankReference = rankReference,
                rankLowDim = rankLowDim,
                kTM = k
            )
            results[[dr]][[paste0("Trustworthiness_k", k)]] <- tm

            if (verbose) message("  Calculating continuity, k=", k)
            ct <- calcContinuityFromRank(
                rankReference = rankReference,
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
            intrs <- colSums((rankReference <= k) * (rankLowDim <= k))
            unin <- colSums(sign((rankReference <= k) + (rankLowDim <= k)))
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
        qt <- coRanking::coranking(
            Xi = as.matrix(distReference),
            X = as.matrix(distLowDim),
            input = "dist"
        )
        lcmc <- coRanking::LCMC(qt)
        ## We follow the implementation in the dimRed package, where Qlocal and
        ## Qglobal are calculated as averages of LCMC values (not Q_NX values)
        Kmax <- which.max(lcmc)
        qlocal <- mean(lcmc[seq(from = 1, to = Kmax, by = 1)])
        qglobal <- mean(lcmc[seq(from = Kmax + 1, to = length(lcmc), by = 1)])
        results[[dr]][["coRankingQlocal"]] <- qlocal
        results[[dr]][["coRankingQglobal"]] <- qglobal
        results[[dr]][["KmaxLCMC"]] <- Kmax
    }

    plots$distcdf <-
        ggplot2::ggplot(
            distdf %>% tidyr::gather(key = "dr", value = "distance"),
            ggplot2::aes(x = distance, color = dr)
        ) +
        ggplot2::stat_ecdf() + ggplot2::theme_bw()

    results <- do.call(dplyr::bind_rows, lapply(names(results), function(nm) {
        data.frame(
            Method = nm, as.data.frame(results[[nm]]),
            stringsAsFactors = FALSE
        )
    }))

    list(scores = results, plots = plots)
}
