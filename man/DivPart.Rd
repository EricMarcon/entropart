\name{DivPart}
\alias{DivPart}
\alias{is.DivPart}
\alias{plot.DivPart}
\alias{autoplot.DivPart}
\alias{summary.DivPart}
\title{
  Diversity Partition of a metacommunity
}
\description{
  Partitions the diversity of a metacommunity into alpha and beta components.
}
\usage{
DivPart(q = 1, MC, Biased = TRUE, Correction = "Best", Tree = NULL,
  Normalize = TRUE, Z = NULL, CheckArguments = TRUE)
is.DivPart(x)
\method{plot}{DivPart}(x, \dots)
\method{autoplot}{DivPart}(object, col = ggplot2::GeomRect$default_aes$fill, 
  border = ggplot2::GeomRect$default_aes$colour, \dots)
\method{summary}{DivPart}(object, \dots)
}
\arguments{
  \item{q}{
  A number: the order of diversity. Default is 1.
  }
  \item{MC}{
  A \code{\link{MetaCommunity}} object.
  }
  \item{Biased}{
  Logical; if \code{FALSE}, a bias correction is appplied.
  }
  \item{Correction}{
  A string containing one of the possible corrections.
	
  The correction must be accepted by \code{\link{AlphaEntropy}}, \code{\link{BetaEntropy}} and \code{\link{GammaEntropy}}. \code{"Best"} is the default value.
  }
  \item{Tree}{
  An object of class \code{\link{hclust}}, "phylo" (see \code{\link[ape]{read.tree}}), \code{\link[ade4]{phylog}} or \code{\link{PPtree}}. The tree must be ultrametric.
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
  \item{x}{
  An object to be tested or plotted.
  }
  \item{object}{
  A \code{MCdiversity} object to be summarized or plotted.
  }
  \item{col}{
  The color used to fill the bars. See "Color Specification" in \code{\link{par}}.
  }
  \item{border}{
  The color of the borders around the bars. See \code{\link{rect}}.
  }
  \item{\dots}{
  Additional arguments to be passed to the generic methods.
  }
}
\details{
  \code{DivPart} partitions the diversity of the metacommunity into alpha and beta components. It supports estimation-bias correction.
  
  If \code{Tree} is provided, the phylogenetic diversity is calculated else if \code{Z} is not \code{NULL}, then similarity-based entropy is calculated.
  
  Beta diversity/entropy is calculated from Gamma and Alpha when bias correction is required, so community values are not available.
}
\value{
  A \code{DivPart} object. It is a list:
  \item{MetaCommunity}{The name of the \code{\link{MetaCommunity}} object containing inventory data.}
  \item{Order}{The value of \code{q}.}
  \item{Biased}{Logical. If \code{FALSE}, bias corrected values of diversity have been computed.}
  \item{Correction}{The estimation bias correction used to calculate diversity.}
  \item{Method}{The method used to calculate entropy ("HCDT", "Similarity-based").}
  \item{Tree}{The phylogenetic or functional tree used to calculate phylodiversity.}
  \item{Normalized}{Logical. Indicates whether phylodiversity is normalized or proportional to the height of the tree.}
  \item{Z}{The matrix used to calculate similarity-based entropy.}
  \item{TotalAlphaDiversity}{The alpha diversity of communities.}
  \item{TotalBetaDiversity}{The beta diversity of communities.}
  \item{GammaDiversity}{The gamma diversity of the metacommunity.}
  \item{CommunityAlphaDiversities}{A vector containing the alpha diversity of each community.}
  \item{TotalAlphaEntropy}{The alpha entropy of communities.}
  \item{TotalBetaEntropy}{The beta entropy of communities.}
  \item{GammaEntropy}{The gamma entropy of the metacommunity.}
  \item{CommunityAlphaEntropies}{A vector containing the alpha entropy of each community.}
  
  \code{DivPart} objects can be summarized and plotted.
}
\references{
  Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The Decomposition of Shannon's Entropy and a Confidence Interval for Beta Diversity. \emph{Oikos} 121(4): 516-522.
  
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.
  
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.

  Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of similarity-based diversity and its bias correction. \emph{HAL} hal-00989454(version 3).
}
\author{
  Eric Marcon <Eric.Marcon@agroparistech.fr>, Bruno Herault <Bruno.Herault@cirad.fr>
}
\seealso{
  \code{\link{DivProfile}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Estimate Shannon diversity.
summary(DivPart(q = 1, Paracou618.MC, Biased = FALSE) -> dp)
plot(dp)
}
