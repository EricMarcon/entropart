\name{Preprocess.Tree}
\alias{Preprocess.Tree}
\title{
  Preprocessing of a phylogenetic tree
}
\description{
  Calculates cuts and intervals of a phylogenetic tree and make it available both in \code{\link{hclust}} and \code{\link{phylo}} formats.
}
\usage{
Preprocess.Tree(Tree)
}
\arguments{
  \item{Tree}{
  An object of class \code{\link{hclust}}, \code{\link{phylo}} or \code{\link[ade4]{phylog}}. The tree must be ultrametric.
  }
}
\details{
  This is an internal function. It is used by \code{\link{PhyloApply}} to obtain values of intervals.
}
\value{
  An object of class \code{\link{PPtree}}.
}
\keyword{internal}
