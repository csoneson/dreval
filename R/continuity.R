#' Calculate continuity based on distance matrices
#'
#' The continuity was proposed by Venna and Kaski, as a local quality measure of
#' a low-dimensional representation. The metric focuses on the preservation of
#' local neighborhoods, and compares the neighborhoods of points in the
#' low-dimensional representation to those in the reference data. Hence, the
#' continuity measure indicates to which degree we can trust that the points
#' closest to a given sample in the reference data set are placed close to the
#' sample also in the low-dimensional representation. The \code{kTM} parameter
#' defines the size of the neighborhoods to consider.
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
#' @return The continuity value.
#'
#' @param distReference N x N matrix or dist object, representing pairwise sample
#'   distances based on the reference (high-dimensional) observed values.
#' @param rankLowDim N x N matrix or dist object, representing pairwise sample
#'   distances based on the low-dimensional representation.
#' @param kTM The number of nearest neighbors
#'
calcContinuityFromDist <- function(distReference, distLowDim, kTM) {
    distReference <- as.matrix(distReference)
    distLowDim <- as.matrix(distLowDim)

    stopifnot(ncol(distReference) == ncol(distLowDim))
    stopifnot(nrow(distReference) == nrow(distLowDim))

    rankReference <- apply(
        as.matrix(distReference), 2, function(w) order(order(w)))
    rankLowDim <- apply(
        as.matrix(distLowDim), 2, function(w) order(order(w)))

    calcContinuityFromRank(
        rankReference = rankReference,
        rankLowDim = rankLowDim, kTM = kTM
    )
}

#' Calculate continuity based on sample rankings
#'
#' The continuity was proposed by Venna and Kaski, as a local quality measure of
#' a low-dimensional representation. The metric focuses on the preservation of
#' local neighborhoods, and compares the neighborhoods of points in the
#' low-dimensional representation to those in the reference data. Hence, the
#' continuity measure indicates to which degree we can trust that the points
#' closest to a given sample in the reference data set are placed close to the
#' sample also in the low-dimensional representation. The \code{kTM} parameter
#' defines the size of the neighborhoods to consider.
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
#' @return The continuity value.
#'
#' @param rankReference N x N matrix, each row/column corresponding to one
#'   sample. The value of entry (i, j) represents the position of sample i in
#'   the ranking of all samples with respect to their distance from sample j,
#'   based on the reference (high-dimensional) observed values. The most similar
#'   sample (i.e., sample j itself) has position 1.
#' @param rankLowDim N x N matrix, each row/column corresponding to one sample.
#'   The value of entry (i, j) represents the position of sample i in the
#'   ranking of all samples with respect to their distance from sample j, based
#'   on the low-dimensional representation. The most similar sample (i.e.,
#'   sample j itself) has position 1.
#' @param kTM The number of nearest neighbors
#'
calcContinuityFromRank <- function(rankReference, rankLowDim, kTM) {
    stopifnot(ncol(rankReference) == ncol(rankLowDim))
    stopifnot(nrow(rankReference) == nrow(rankLowDim))

    N <- ncol(rankReference)

    1 - 2/(N * kTM * (2 * N - 3 * kTM - 1)) *
        sum(vapply(seq_len(ncol(rankReference)), function(i) {
            sum((rankLowDim[, i] - kTM) * (rankReference[, i] <= kTM) *
                    (rankLowDim[, i] > kTM))
        }, NA_real_))
}
