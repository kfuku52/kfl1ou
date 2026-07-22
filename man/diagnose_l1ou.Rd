\name{diagnose_l1ou}
\alias{diagnose_l1ou}
\title{Diagnose numerical and residual problems in a fitted model}
\usage{diagnose_l1ou(model)}
\arguments{\item{model}{Fitted object of class \code{"l1ou"}.}}
\value{A list containing optimizer, covariance, standardized-residual, and
potential tip-outlier diagnostics.}
\description{Checks convergence and boundaries, covariance conditioning,
whitened residuals, residual cross-trait correlation, and unusually large
tip residuals.}
