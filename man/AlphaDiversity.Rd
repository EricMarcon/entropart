\name{AlphaDiversity}
\alias{AlphaDiversity}
\title{
  Reduced-bias alpha diversity of a metacommunity
}
\description{
  Calculates the eeduced-bias total alpha diversity of order \eqn{q} of communities.
}
\usage{
AlphaDiversity(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, 
  Z = NULL, CheckArguments = TRUE)
}
\arguments{
  \item{MC}{
    A \code{\link{MetaCommunity}} object.
  }
  \item{q}{
    A number: the order of diversity. Default is 1 for Shannon diversity.
  }
  \item{Correction}{
    A string containing one of the possible corrections accepted by \code{\link{AlphaEntropy}} or \code{"None"} or \code{"Best"}, the default value.
  }
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{Normalize}{
    If \code{TRUE} (default), diversity is not affected by the height of the tree.\cr
    If \code{FALSE}, diversity is proportional to the height of the tree.
  }
  \item{Z}{
    A relatedness matrix, \emph{i.e.} a square matrix whose terms are all positive, strictly positive on the diagonal. Generally, the matrix is a similarity matrix, \emph{i.e.} the diagonal terms equal 1 and other terms are between 0 and 1.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Entropy is calculated by \code{\link{AlphaEntropy}} and transformed into diversity.
}
\value{
  An \code{\link{MCdiversity}} object containing diversity values of each community and of the metacommunity.
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
  \code{\link{AlphaEntropy}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Calculate Simpson alpha diversity
  summary(AlphaDiversity(Paracou618.MC, 2))
  # Compare without correction
  summary(AlphaDiversity(Paracou618.MC, 2, Correction = "None"))
  # Estimate phylogenetic Simpson alpha diversity
  summary(AlphaDiversity(Paracou618.MC, 2, Tree = Paracou618.Taxonomy) -> e)
  plot(e)  
}
