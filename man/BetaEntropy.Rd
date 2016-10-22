\name{BetaEntropy}
\alias{BetaEntropy}
\title{
  Reduced-bias beta entropy of a metacommunity
}
\description{
  Calculates the reduced-bias beta entropy of order \eqn{q} between communities.
}
\usage{
BetaEntropy(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, 
  Z = NULL, CheckArguments = TRUE)
}
\arguments{
  \item{MC}{
    A \code{\link{MetaCommunity}} object.
  }
  \item{q}{
    A number: the order of diversity. Default is 1 for Shannon entropy.
  }
  \item{Correction}{
    A string containing one of the possible corrections accepted by the bias-corrected entropy function (see details) or \code{"None"} or \code{"Best"}, the default value.
  }
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{Normalize}{
    If \code{TRUE} (default), the entropy returned by the function is normalized by the height of the tree (it is the weighted average value of the entropy in each slice).\cr
    If \code{FALSE}, it is the unnormalized weighted sum of the results.
  }
  \item{Z}{
    A relatedness matrix, \emph{i.e.} a square matrix whose terms are all positive, strictly positive on the diagonal. Generally, the matrix is a similarity matrix, \emph{i.e.} the diagonal terms equal 1 and other terms are between 0 and 1.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  If \code{Tree} is not \code{NULL}, then phylogenetic entropy is calculated by \code{\link{bcPhyloBetaEntropy}}; else, if \code{Z} is not \code{NULL}, then similarity-based entropy is calculated by \code{\link{bcHqzBeta}}; else, neutral entropy is calculated by \code{\link{bcTsallisBeta}}.
  
  The reduced-bias beta entropy of each community is calculated and summed according to community weights.

  Note that beta entropy is related to alpha entropy (if \eqn{q} is not 1) and cannot be compared accross communities (Jost, 2007). Do rather calculate the \code{\link{BetaDiversity}} of the metacommunity.
}
\value{
  An \code{\link{MCentropy}} object containing entropy values of each community and of the metacommunity.
}
\references{
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.
  
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.

  Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of similarity-based diversity and its bias correction. \emph{HAL} hal-00989454(version 1).
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{bcTsallisBeta}}, \code{\link{BetaDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Estimate Shannon beta entropy
  summary(BetaEntropy(Paracou618.MC, 1))
  # Compare without correction
  summary(BetaEntropy(Paracou618.MC, 1, Correction = "None"))
  # Estimate phylogenetic Shannon beta entropy
  summary(BetaEntropy(Paracou618.MC, 1, Tree = Paracou618.Taxonomy) -> e)
  plot(e)  
}
