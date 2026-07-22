\name{confint.l1ou}
\alias{confint.l1ou}
\title{Confidence intervals for fitted kfl1ou parameters}
\usage{\method{confint}{l1ou}(object, parm, level = 0.95,
  method = c("parametric", "wald"), nsim = 100L, seed = NULL, ...)}
\arguments{
\item{object}{Fitted object of class \code{"l1ou"}.}
\item{parm}{Parameter names or indices; all reported parameters by default.}
\item{level}{Confidence level.}
\item{method}{Parametric-bootstrap or Hessian-based Wald intervals.}
\item{nsim}{Number of parametric bootstrap replicates.}
\item{seed}{Optional random seed.}
\item{...}{Additional arguments, currently ignored.}
}
\value{A two-column matrix of confidence limits. Failed parametric-bootstrap
refits are omitted and counted in attributes.}
\description{Quantifies parameter uncertainty by fixed-configuration
parametric refitting. Wald intervals are available for alpha when a suitable
optimizer Hessian is available.}
