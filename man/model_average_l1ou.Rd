\name{model_average_l1ou}
\alias{model_average_l1ou}
\title{Information-criterion model averaging over shift configurations}
\usage{model_average_l1ou(model, delta.max = 10, max.models = Inf)}
\arguments{
\item{model}{Fitted \code{"l1ou"} object with a search profile.}
\item{delta.max}{Largest score difference retained from the best model.}
\item{max.models}{Optional cap on retained models.}
}
\value{Normalized weights, refitted models, averaged means, optima and process
parameters, and edge-inclusion weights.}
\description{Propagates shift-configuration uncertainty into predictions and
parameter summaries using normalized information-criterion weights.}
