#' Calculate continuity based on distance matrices
#'
#' The continuity was proposed by Kaski et al, as a local quality measure
#' of a low-dimensional representation. The metric focuses on the preservation
#' of local neighborhoods, and compares the neighborhoods of points in the
#' low-dimensional representation to those in the original data. Hence, the
#' continuity measure indicates to which degree we can trust that the points
#' closest to a given sample in the original data set are placed close to the
#' sample also in the low-dimensional representation. The \code{kTM} parameter
#' defines the size of the neighborhoods to consider.
#'
#' @references Kaski, S., Nikkilä, J., Oja, M., Venna, J., Törönen, P., and
#'   Castrén, E. (2003). Trustworthiness and metrics in visualizing similarity
#'   of gene expression. BMC Bioinformatics 4:48.
#'
#' @keywords internal
#'
#' @author Charlotte Soneson
#'
#' @return The continuity value.
#'
#' @param distOriginal N x N matrix or dist object, representing pairwise sample
#'   distances based on the original (high-dimensional) observed values.
#' @param rankLowDim N x N matrix or dist object, representing pairwise sample
#'   distances based on the low-dimensional representation.
#' @param kTM The number of nearest neighbors
#'
calcContinuityFromDist <- function(distOriginal, distLowDim, kTM) {
    distOriginal <- as.matrix(distOriginal)
    distLowDim <- as.matrix(distLowDim)

    stopifnot(ncol(distOriginal) == ncol(distLowDim))
    stopifnot(nrow(distOriginal) == nrow(distLowDim))

    rankOriginal <- apply(
        as.matrix(distOriginal), 2, function(w) order(order(w)))
    rankLowDim <- apply(
        as.matrix(distLowDim), 2, function(w) order(order(w)))

    calcContinuityFromRank(
        rankOriginal = rankOriginal,
        rankLowDim = rankLowDim, kTM = kTM
    )
}

#' Calculate continuity based on sample rankings
#'
#' The continuity was proposed by Kaski et al, as a local quality measure
#' of a low-dimensional representation. The metric focuses on the preservation
#' of local neighborhoods, and compares the neighborhoods of points in the
#' low-dimensional representation to those in the original data. Hence, the
#' continuity measure indicates to which degree we can trust that the points
#' closest to a given sample in the original data set are placed close to the
#' sample also in the low-dimensional representation. The \code{kTM} parameter
#' defines the size of the neighborhoods to consider.
#'
#' @references Kaski, S., Nikkilä, J., Oja, M., Venna, J., Törönen, P., and
#'   Castrén, E. (2003). Trustworthiness and metrics in visualizing similarity
#'   of gene expression. BMC Bioinformatics 4:48.
#'
#' @keywords internal
#'
#' @author Charlotte Soneson
#'
#' @return The continuity value.
#'
#' @param rankOriginal N x N matrix, each row/column corresponding to one
#'   sample. The value of entry (i, j) represents the position of sample i in
#'   the ranking of all samples with respect to their distance from sample j,
#'   based on the original (high-dimensional) observed values. The most similar
#'   sample (i.e., sample j itself) has position 1.
#' @param rankLowDim N x N matrix, each row/column corresponding to one sample.
#'   The value of entry (i, j) represents the position of sample i in the
#'   ranking of all samples with respect to their distance from sample j, based
#'   on the low-dimensional representation. The most similar sample (i.e.,
#'   sample j itself) has position 1.
#' @param kTM The number of nearest neighbors
#'
calcContinuityFromRank <- function(rankOriginal, rankLowDim, kTM) {
    stopifnot(ncol(rankOriginal) == ncol(rankLowDim))
    stopifnot(nrow(rankOriginal) == nrow(rankLowDim))

    N <- ncol(rankOriginal)

    1 - 2/(N * kTM * (2 * N - 3 * kTM - 1)) *
        sum(vapply(seq_len(ncol(rankOriginal)), function(i) {
            sum((rankLowDim[, i] - kTM) * (rankOriginal[, i] <= kTM) *
                    (rankLowDim[, i] > kTM))
        }, NA_real_))
}
