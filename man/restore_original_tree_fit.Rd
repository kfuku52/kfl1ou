% Generated manually: do not edit by hand if roxygen is introduced later
\name{restore_original_tree_fit}
\alias{restore_original_tree_fit}
\title{Restore a pruned fit to the original tree}
\usage{
restore_original_tree_fit(
  model,
  original.tree,
  representative = c("tipward", "rootward")
)
}
\arguments{
\item{model}{object of class \code{"l1ou"} returned by \pkg{kfl1ou}.}

\item{original.tree}{original tree before pruning tips with missing data.}

\item{representative}{how to choose a single representative original edge when
a pruned edge maps to a collapsed path in \code{original.tree}. The default
\code{"tipward"} picks the deepest edge on that path; \code{"rootward"}
picks the shallowest one.}
}
\value{
An object of class \code{"restored_l1ou"} and \code{"l1ou"} with the
original tree, tip-level matrices padded back to the original tip set, and a
\code{restoration} component containing:

\item{removed.tips}{tip labels that were absent from the fitted tree.}
\item{pruned.shift.configuration}{shift edges on the fitted tree.}
\item{shift.edge.paths}{list mapping each fitted shift edge to the
corresponding path of edge indices in \code{original.tree}.}
\item{original.to.pruned.edge}{integer vector mapping each edge of
\code{original.tree} back to the corresponding fitted-tree edge, or \code{NA}
for edges supported only by removed tips.}
\item{pruned.to.original.edge}{list mapping each fitted-tree edge to the set of
original-tree edges it represents.}
}
\description{
This helper is intended for reporting after a fit was run on a tree pruned for
missing-tip handling. It pads tip-level outputs back to the original tip set
and maps each shift edge to the corresponding edge path in the original tree.
}
\details{
When pruning drops tips, some edges in the fitted tree correspond to collapsed
paths in the original tree. In those cases, \code{shift.configuration} in the
restored object uses representative edges chosen according to
\code{representative}, while the full ambiguous paths are retained in
\code{$restoration$shift.edge.paths}.
}
\examples{
data(lizard.tree, lizard.traits)
keep <- lizard.tree$tip.label[1:8]
tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
Y <- as.matrix(lizard.traits[keep, 1, drop = FALSE])
Y[1, 1] <- NA_real_

adjusted <- adjust_data(tree, Y, normalize = FALSE)
fit <- estimate_shift_configuration(adjusted$tree, adjusted$Y, max.nShifts = 1)
restored <- restore_original_tree_fit(fit, tree)
print(restored)
}
