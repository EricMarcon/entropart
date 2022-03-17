\name{ChaoPD}
\alias{ChaoPD}
\title{
  Phylogenetic Diversity of a Community
}
\description{
  Calculates the phylogenetic diversity of order \eqn{q} of a probability vector.
}
\usage{
ChaoPD(Ps, q = 1, PhyloTree, Normalize = TRUE, Prune = FALSE, CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
  A probability vector, summing to 1.
  }
  \item{q}{
  A number: the order of diversity. Default is 1.
  }
  \item{PhyloTree}{
  An object of class \code{\link{hclust}}, \code{\link{phylo}}, \code{\link[ade4]{phylog}} or \code{\link{PPtree}}. The tree is not necessarily ultrametric.
  }
  \item{Normalize}{
  If \code{TRUE} (default), diversity is not affected by the height of the tree.\cr
  If \code{FALSE}, it is proportional to the height of the tree.
  }
  \item{Prune}{
  What to do when somes species are in the tree but not in \code{Ps}?\cr
  If \code{TRUE}, the tree is pruned to keep species of \code{Ps} only. The height of the tree may be changed if a pruned branch is related to the root.\cr
  If \code{FALSE} (default), species with probability 0 are added in \code{Ps}.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The phylogenetic diversity is calculated following Chao \emph{et al.} (2010). The result is identical to the total diversity calculated by \code{\link{PhyloDiversity}} but it is much faster. A single value is returned instead of a \code{\link{PhyloDiversity}} object, and no bias correction is available.
  
  The \code{Normalize} arguments allows calculating either \eqn{^{q}\bar{D}(T)} (if \code{TRUE}) or \eqn{^{q}PD(T)} if \code{FALSE}.

  Diversity can be calculated for non ultrametric trees following Chao \emph{et al.} (2010) even though the meaning of the result is not so clear (Leinster and Cobold, 2011).
}
\value{
  A named number equal the diversity of the community. The name is "None" to recall that no bias correction is available.
}
\references{
  Chao, A., Chiu, C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B} 365(1558): 3599-609.
  
  Leinster, T. and Cobbold, C. (2011). Measuring diversity: the importance of species similarity. \emph{Ecology} 93(3): 477-489.
}
\seealso{
  \code{\link{PhyloDiversity}}, \code{\link{AllenH}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- Paracou618.MC$Ps

# Calculate the phylogenetic Simpson diversity of the plot
(ChaoPD(Paracou618.MC$Ps, 2, Paracou618.Taxonomy, Normalize=TRUE))

# Calculate it using PhyloDiversity
# (more powerful but much slower if the tree has many periods)
PhyloDiversity(Paracou618.MC$Ps, 2, Paracou618.Taxonomy, Normalize=TRUE) -> phyD
summary(phyD)
}
