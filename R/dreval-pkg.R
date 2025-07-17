#' dreval
#'
#' dreval evaluates and compares multiple reduced dimension
#' representations, based on how well they retain structure from
#' the original data set.
#'
#' @keywords internal
"_PACKAGE"

globalVariables(c(
    "reference", "lowdim", "stat", "density",
    "distance", "dimensionality", "KSStatDist",
    "EuclDistBetweenDists", "SammonStress",
    "Method", "metric", "score", "total"))
