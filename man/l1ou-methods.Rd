\name{l1ou-methods}
\alias{l1ou-methods}
\alias{coef.l1ou}
\alias{fitted.l1ou}
\alias{residuals.l1ou}
\alias{logLik.l1ou}
\alias{nobs.l1ou}
\alias{vcov.l1ou}
\alias{predict.l1ou}
\title{Standard methods for fitted kfl1ou models}
\usage{
\method{coef}{l1ou}(object, ...)
\method{fitted}{l1ou}(object, ...)
\method{residuals}{l1ou}(object, ...)
\method{logLik}{l1ou}(object, ...)
\method{nobs}{l1ou}(object, ...)
\method{vcov}{l1ou}(object, ...)
\method{predict}{l1ou}(object, type = c("response", "optimum"), ...)
}
\arguments{
\item{object}{Fitted object of class \code{"l1ou"}.}
\item{type}{Prediction type: expected tip response or adaptive optimum.}
\item{...}{Additional arguments, currently ignored.}
}
\value{The conventional extracted model component. \code{vcov} returns the
fitted evolutionary innovation covariance among traits.}
\description{Extract coefficients, fitted values, residuals, likelihood,
sample size, evolutionary covariance, or fitted tip predictions.}
