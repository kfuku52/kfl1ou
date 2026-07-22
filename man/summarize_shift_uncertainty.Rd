\name{summarize_shift_uncertainty}
\alias{summarize_shift_uncertainty}
\title{Summarize uncertainty in bootstrap shift locations}
\usage{summarize_shift_uncertainty(bootstrap.result, tree = NULL)}
\arguments{
\item{bootstrap.result}{Result from \code{\link{l1ou_bootstrap_support}}.}
\item{tree}{Optional phylogeny used to label edges.}
}
\value{Configuration frequencies, edge inclusion probabilities, edge
co-selection probabilities, and the distribution of shift counts.}
\description{Turns bootstrap shift configurations into stability summaries
that retain both marginal edge support and dependence among selected edges.}
