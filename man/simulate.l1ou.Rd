\name{simulate.l1ou}
\alias{simulate.l1ou}
\title{Simulate from a fitted kfl1ou model}
\usage{\method{simulate}{l1ou}(object, nsim = 1L, seed = NULL,
  preserve.missing = TRUE, ...)}
\arguments{
\item{object}{Fitted object of class \code{"l1ou"}.}
\item{nsim}{Number of simulated data sets.}
\item{seed}{Optional random seed.}
\item{preserve.missing}{Retain the fitted data's missingness pattern.}
\item{...}{Additional arguments, currently ignored.}
}
\value{A named list of simulated trait matrices.}
\description{Generates Gaussian OU data using the fitted phylogenetic,
cross-trait, and observation-error covariance.}
