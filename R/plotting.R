#' Plot summary of the dimension reduction method ranking across metrics
#'
#' For each metric, rank the evaluated reduced dimension representations by
#' performance, and plot a summary of the overall ranking. Metrics evaluating
#' local and global structure preservations are colored in red and blue,
#' respectively.
#'
#' @param dreSummary A \code{data.frame} with the values of the evaluation
#'   metrics, typically the \code{"scores"} element of the output of
#'   \code{dreval()}.
#' @param metrics A character vector with the metrics to include in the summary.
#'   Must be a subset of the column names of \code{dreSummary}. If \code{NULL},
#'   all metrics will be used.
#'
#' @author Charlotte Soneson
#'
#' @return Nothing is returned, but a plot is generated.
#'
#' @export
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr select mutate vars
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw scale_fill_manual
#'
plotRankSummary <- function(dreSummary, metrics = NULL) {
    myfun <- function(w) {
        order(order(w))
    }

    ## Define metrics to use
    cn <- colnames(dreSummary)
    if (is.null(metrics)) {
        metrics <- cn
    }
    global <- intersect(
        metrics,
        c("coRankingQglobal", "EuclDistBetweenDists", "KSStatDist",
          "PearsonCorrDist", "SammonStress", "SpearmanCorrDist")
    )
    local <- intersect(
        metrics,
        c(grep("Continuity_k", cn, value = TRUE),
          grep("Trustworthiness_k", cn, value = TRUE),
          "coRankingQlocal",
          grep("MeanJaccard_k", cn, value = TRUE))
    )

    ## Define colors
    globalcols <- grDevices::colorRampPalette(
        c("darkblue", "lightblue"))(length(global))
    names(globalcols) <- global
    localcols <- grDevices::colorRampPalette(
        c("darkred", "pink"))(length(local))
    names(localcols) <- local

    dreSummary <- dreSummary %>% dplyr::select(-dimensionality) %>%
        dplyr::mutate(KSStatDist = -KSStatDist) %>%
        dplyr::mutate(EuclDistBetweenDists = -EuclDistBetweenDists) %>%
        dplyr::mutate(SammonStress = -SammonStress) %>%
        dplyr::select(c("Method", global, local)) %>%
        dplyr::mutate_at(dplyr::vars(-Method), myfun) %>%
        tidyr::gather(key = "metric", value = "score", -Method) %>%
        dplyr::mutate(metric = factor(metric, levels = c(global, local)))

    ggplot2::ggplot(dreSummary,
                    ggplot2::aes(x = Method, y = score, fill = metric)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = c(globalcols, localcols)) +
        ggplot2::labs(y = "Total rank score (high = good)")
}
