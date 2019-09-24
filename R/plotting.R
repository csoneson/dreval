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
#'   all metrics will be used. It can also be "global" or "local", in which case
#'   all the global or local metrics, respectively, will be used.
#' @param sortBars A character scalar indicating whether/how to sort the bars in
#'   the output. Either "decreasing", "increasing" or "none" (in which case the
#'   input order will be used).
#' @param scoreType A character scalar indicating what type of values to show in
#'   the plot. Either "rank" or "rescale". If set to "rank", the representations
#'   will be ranked for each metric (with the best one assigned the highest
#'   rank). If set to "rescale", the scores for each metric will first, if
#'   necessary, be inverted so that a high (positive) value corresponds to
#'   better performance, and then be linearly rescaled, mapping the lowest score
#'   to 1 and the highest to P, where P is the number of evaluated
#'   representations. If the original scores are approximately equally spaced
#'   between the highest and lowest observed values, this gives similar results
#'   as setting \code{scoreType} to "rank". However, if some of the scores are
#'   very similar to each other, the "rescale" approach allows them to get a
#'   similar rank score rather than forcing a uniform difference between
#'   successive scores.
#'
#' @author Charlotte Soneson
#'
#' @return Nothing is returned, but a plot is generated.
#'
#' @export
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr select mutate vars group_by ungroup arrange mutate_at desc
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw scale_fill_manual
#'
plotRankSummary <- function(dreSummary, metrics = NULL,
                            sortBars = "decreasing", scoreType = "rank") {
    scorefun <- function(w, scoreType = "rank") {
        if (scoreType == "rank") {
            return(order(order(w)))
        } else if (scoreType == "rescale") {
            m <- min(w)
            M <- max(w)
            P <- length(w)
            return(((P - 1) * w + M - m * P)/(M - m))
        } else {
            stop("Unknown 'scoreType'")
        }
    }

    ## Define metrics to use
    cn <- colnames(dreSummary)

    global <- c("coRankingQglobal", "EuclDistBetweenDists", "KSStatDist",
                "PearsonCorrDist", "SammonStress", "SpearmanCorrDist")
    local <- c(grep("Continuity_k", cn, value = TRUE),
               grep("Trustworthiness_k", cn, value = TRUE),
               "coRankingQlocal",
               grep("MeanJaccard_k", cn, value = TRUE))

    if (is.null(metrics)) {
        metrics <- cn
    } else if (length(metrics) == 1 && metrics == "global") {
        metrics <- intersect(cn, global)
    } else if (length(metrics) == 1 && metrics == "local") {
        metrics <- intersect(cn, local)
    }

    global <- intersect(metrics, global)
    local <- intersect(metrics, local)

    ## Define colors
    globalcols <- grDevices::colorRampPalette(
        c("darkblue", "lightblue"))(length(global))
    names(globalcols) <- global
    localcols <- grDevices::colorRampPalette(
        c("darkred", "pink"))(length(local))
    names(localcols) <- local

    ## Make sure that for all scores, high values represent good performance
    dreSummary <- dreSummary %>% dplyr::select(-dimensionality) %>%
        dplyr::mutate(KSStatDist = -KSStatDist) %>%
        dplyr::mutate(EuclDistBetweenDists = -EuclDistBetweenDists) %>%
        dplyr::mutate(SammonStress = -SammonStress) %>%
        dplyr::select(c("Method", global, local)) %>%
        dplyr::mutate_at(dplyr::vars(-Method), scorefun, scoreType = scoreType) %>%
        tidyr::gather(key = "metric", value = "score", -Method) %>%
        dplyr::mutate(metric = factor(metric, levels = c(global, local))) %>%
        dplyr::ungroup()

    if (sortBars == "decreasing") {
        dreSummary <- dreSummary %>% dplyr::group_by(Method) %>%
            dplyr::mutate(total = sum(score)) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(dplyr::desc(total)) %>%
            dplyr::mutate(Method = factor(Method, levels = unique(Method))) %>%
            dplyr::select(-total)
    } else if (sortBars == "increasing") {
        dreSummary <- dreSummary %>% dplyr::group_by(Method) %>%
            dplyr::mutate(total = sum(score)) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(total) %>%
            dplyr::mutate(Method = factor(Method, levels = unique(Method))) %>%
            dplyr::select(-total)
    } else if (sortBars == "none") {
        dreSummary <- dreSummary %>%
            dplyr::mutate(Method = factor(Method, levels = unique(Method)))
    } else {
        stop("'sortBars' must be either 'decreasing', 'increasing' or 'none'.")
    }

    ggplot2::ggplot(dreSummary,
                    ggplot2::aes(x = Method, y = score, fill = metric)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = c(globalcols, localcols)) +
        ggplot2::labs(y = "Total rank score (high = good)")
}
