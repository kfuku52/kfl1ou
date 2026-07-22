#' kfl1ou: OU shift detection and comparative trait modelling
#'
#' Internal package declarations used to generate the namespace.
#'
#' @keywords internal
#' @importFrom ape Nedge drop.tip edgelabels is.binary is.ultrametric multi2di plot.phylo
#' @importFrom grDevices rainbow
#' @importFrom graphics axis barplot layout mtext par
#' @importFrom stats AIC coef confint delete.response fitted formula logLik lm.fit
#' @importFrom stats model.frame model.matrix model.response nobs optim optimHess
#' @importFrom stats optimize printCoefmat predict profile pt quantile reorder
#' @importFrom stats residuals sd setNames simulate terms var vcov
#' @importFrom utils capture.output getFromNamespace
#' @importFrom Rcpp evalCpp
#' @useDynLib kfl1ou, .registration = TRUE
"_PACKAGE"
