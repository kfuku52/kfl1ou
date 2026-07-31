#' kfl1ou: OU shift detection and comparative trait modelling
#'
#' Fits Ornstein-Uhlenbeck comparative models to detect changes in adaptive
#' optima on a phylogeny. The package supports univariate and multivariate
#' traits, convergent regimes, measurement error, full evolutionary trait
#' covariance, bootstrap uncertainty, model averaging, and tree sensitivity.
#'
#' @section Getting started:
#' Use \code{\link[=adjust_data]{adjust_data()}} to validate and align a tree
#' and trait matrix, then call
#' \code{\link[=estimate_shift_configuration]{estimate_shift_configuration()}}
#' to search for shifts or \code{\link[=fit_OU]{fit_OU()}} to fit a specified
#' configuration. Use \code{\link[=diagnose_l1ou]{diagnose_l1ou()}},
#' \code{\link[=l1ou_bootstrap_support]{l1ou_bootstrap_support()}}, and
#' \code{\link[=sensitivity_l1ou]{sensitivity_l1ou()}} before interpreting a
#' selected model. See
#' \code{vignette("multivariate-inference", package = "kfl1ou")} for a correlated
#' multivariate workflow.
#'
#' @section Model scope and interpretation:
#' Shift locations identify phylogenetically consistent changes in expected
#' trait optima; extant-tip data generally do not identify exact historical
#' shift times. OU adaptation rates, optima, measurement error, and covariance
#' structure can be weakly identified. Report optimizer and boundary
#' diagnostics, the candidate configurations evaluated, and sensitivity to
#' alpha bounds, root models, phylogenetic trees, and resampling choices.
#'
#' @section Citation:
#' Run \code{citation("kfl1ou")} for the package citation and the original l1ou
#' methodology reference.
#'
#' @keywords package
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
