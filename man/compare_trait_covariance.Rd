\name{compare_trait_covariance}
\alias{compare_trait_covariance}
\alias{print.l1ou_covariance_comparison}
\title{Compare diagonal and correlated-trait OU covariance models}
\usage{
compare_trait_covariance(tree, Y, shift.configuration = integer(),
  nboot = 100L, seed = NULL,
  root.model = c("OUfixedRoot", "OUrandomRoot"),
  alpha.structure = c("diagonal", "shared"),
  likelihood.engine = c("auto", "dense", "pruning"),
  optimizer.starts = 5L, measurement_error = FALSE,
  input_error = NULL, quietly = TRUE)
\method{print}{l1ou_covariance_comparison}(x, ...)
}
\arguments{
\item{tree}{Ultrametric phylogeny in postorder.}
\item{Y}{Numeric trait matrix aligned to the tree tips.}
\item{shift.configuration}{Fixed shift edges used in both models.}
\item{nboot}{Number of null simulations; zero skips bootstrap calibration.}
\item{seed}{Optional random seed.}
\item{root.model}{Ancestral-root model.}
\item{alpha.structure}{Full-model alpha structure.}
\item{likelihood.engine}{Full-model likelihood engine.}
\item{optimizer.starts}{Number of optimizer starts.}
\item{measurement_error}{Estimate additional observation variance.}
\item{input_error}{Optional known observation-error variances.}
\item{quietly}{Suppress progress messages.}
\item{x}{Comparison object.}
\item{...}{Additional print arguments.}
}
\value{An object containing both fits, the observed likelihood ratio,
bootstrap statistics, and a calibrated p-value.}
\description{Tests whether cross-trait evolutionary covariances improve fit,
using a parametric-bootstrap likelihood-ratio reference distribution.}
