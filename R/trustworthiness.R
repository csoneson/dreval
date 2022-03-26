#' Calculate trustworthiness based on distance matrices
#'
#' The trustworthiness was proposed by Venna and Kaski, as a local quality
#' measure of a low-dimensional representation. The metric focuses on the
#' preservation of local neighborhoods, and compares the neighborhoods of points
#' in the low-dimensional representation to those in the reference data. Hence,
#' the trustworthiness measure indicates to which degree we can trust that the
#' points placed closest to a given sample in the low-dimensional representation
#' are really close to the sample also in the reference data set. The \code{kTM}
#' parameter defines the size of the neighborhoods to consider.
#'
#' @references
#'   Venna J., Kaski S. (2001). Neighborhood preservation in nonlinear
#'   projection methods: An experimental study. In Dorffner G., Bischof H.,
#'   Hornik K., editors, Proceedings of ICANN 2001, pp 485–491. Springer,
#'   Berlin.
#'
#' @keywords internal
#'
#' @author Charlotte Soneson
#'
#' @return The trustworthiness value.
#'
#' @param distReference N x N matrix or dist object, representing pairwise
#'   sample distances based on the reference (high-dimensional) observed values.
#'   For each column, samples (rows) will be ranked by the provided distances.
#' @param rankLowDim N x N matrix or dist object, representing pairwise sample
#'   distances based on the low-dimensional representation.
#'   For each column, samples (rows) will be ranked by the provided distances.
#' @param kTM The number of nearest neighbors (excluding the sample itself).
#'
calcTrustworthinessFromDist <- function(distReference, distLowDim, kTM) {
    distReference <- as.matrix(distReference)
    distLowDim <- as.matrix(distLowDim)

    ## Check that matrices have the same dimensions, and that the
    ## observations are ordered in the same way
    stopifnot(ncol(distReference) == ncol(distLowDim))
    stopifnot(nrow(distReference) == nrow(distLowDim))
    stopifnot(rownames(distReference) == colnames(distReference))
    stopifnot(rownames(distReference) == rownames(distLowDim))
    stopifnot(rownames(distLowDim) == colnames(distLowDim))

    diag(distReference) <- NA
    diag(distLowDim) <- NA

    rankReference <- apply(distReference, 2, function(w) order(order(w)))
    rankLowDim <- apply(distLowDim, 2, function(w) order(order(w)))

    diag(rankReference) <- 0
    diag(rankLowDim) <- 0

    rownames(rankReference) <- rownames(distReference)
    rownames(rankLowDim) <- rownames(distLowDim)

    calcTrustworthinessFromRank(rankReference = rankReference,
                                rankLowDim = rankLowDim, kTM = kTM)
}

#' Calculate trustworthiness based on sample rankings
#'
#' The trustworthiness was proposed by Venna and Kaski, as a local quality
#' measure of a low-dimensional representation. The metric focuses on the
#' preservation of local neighborhoods, and compares the neighborhoods of points
#' in the low-dimensional representation to those in the reference data. Hence,
#' the trustworthiness measure indicates to which degree we can trust that the
#' points placed closest to a given sample in the low-dimensional representation
#' are really close to the sample also in the reference data set. The \code{kTM}
#' parameter defines the size of the neighborhoods to consider.
#'
#' @references
#'   Venna J., Kaski S. (2001). Neighborhood preservation in nonlinear
#'   projection methods: An experimental study. In Dorffner G., Bischof H.,
#'   Hornik K., editors, Proceedings of ICANN 2001, pp 485–491. Springer,
#'   Berlin.
#'
#' @keywords internal
#'
#' @author Charlotte Soneson
#'
#' @return The trustworthiness value.
#'
#' @param rankReference N x N matrix, each row/column corresponding to one
#'   sample. The value of entry (i, j) represents the position of sample i in
#'   the ranking of all samples with respect to their distance from sample j,
#'   based on the reference (high-dimensional) observed values. The sample
#'   itself has rank 0.
#' @param rankLowDim N x N matrix, each row/column corresponding to one sample.
#'   The value of entry (i, j) represents the position of sample i in the
#'   ranking of all samples with respect to their distance from sample j, based
#'   on the low-dimensional representation. The sample
#'   itself has rank 0.
#' @param kTM The number of nearest neighbors.
#'
calcTrustworthinessFromRank <- function(rankReference, rankLowDim, kTM) {
    stopifnot(ncol(rankReference) == ncol(rankLowDim))
    stopifnot(nrow(rankReference) == nrow(rankLowDim))
    stopifnot(rownames(rankReference) == colnames(rankReference))
    stopifnot(rownames(rankReference) == rownames(rankLowDim))
    stopifnot(rownames(rankLowDim) == colnames(rankLowDim))

    N <- ncol(rankReference)

    if (kTM < N/2) {
        normConst <- N * kTM * (2 * N - 3 * kTM - 1)
    } else {
        normConst <- N * (N - kTM) * (N - kTM - 1)
    }

    1 - 2/(normConst) *
        sum(vapply(seq_len(ncol(rankLowDim)), function(i) {
            sum((rankReference[, i] - kTM) * (rankLowDim[, i] <= kTM) *
                    (rankReference[, i] > kTM))
        }, NA_real_))
}
