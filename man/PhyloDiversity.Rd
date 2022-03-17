\name{PhyloDiversity}
\alias{PhyloDiversity}
\alias{bcPhyloDiversity}
\alias{PhyloDiversity.ProbaVector}
\alias{PhyloDiversity.AbdVector}
\alias{PhyloDiversity.integer}
\alias{PhyloDiversity.numeric}
\alias{is.PhyloDiversity}
\alias{summary.PhyloDiversity}
\title{
  Phylogenetic Diversity of a Community
}
\description{
  Calculates the phylogenetic diversity of order \eqn{q} of a probability vector.
}
\usage{
PhyloDiversity(NorP, q = 1, Tree, Normalize = TRUE, \dots)
bcPhyloDiversity(Ns, q = 1, Tree, Normalize = TRUE, Correction = "Best",
  CheckArguments = TRUE)
\method{PhyloDiversity}{ProbaVector}(NorP, q = 1, Tree, Normalize = TRUE,  
  \dots, CheckArguments = TRUE, Ps = NULL)
\method{PhyloDiversity}{AbdVector}(NorP, q = 1, Tree, Normalize = TRUE, 
  Correction = "Best", \dots, CheckArguments = TRUE, Ns = NULL)
\method{PhyloDiversity}{integer}(NorP, q = 1, Tree, Normalize = TRUE, 
  Correction = "Best", \dots, CheckArguments = TRUE, Ns = NULL)
\method{PhyloDiversity}{numeric}(NorP, q = 1, Tree, Normalize = TRUE, 
  Correction = "Best", \dots, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
is.PhyloDiversity(x)
\method{summary}{PhyloDiversity}(object, \dots)
}
\arguments{
  \item{Ps}{
  A probability vector, summing to 1.
  }
  \item{Ns}{
  A numeric vector containing species abundances.
  }
  \item{NorP}{
  A numeric vector, an integer vector, an abundance vector (\code{\link{AbdVector}}) or a probability vector (\code{\link{ProbaVector}}). Contains either abundances or probabilities.
  }
  \item{q}{
  A number: the order of diversity. Default is 1.
  }
  \item{Tree}{
  An object of class \code{\link{hclust}}, \code{\link{phylo}}, \code{\link[ade4]{phylog}} or \code{\link{PPtree}}. The tree must be ultrametric.
  }
  \item{Normalize}{
  If \code{TRUE} (default), the \code{Total} diversity is not affected by the height of the tree.\cr
  If \code{FALSE}, it is proportional to the height of the tree.\cr
  }
  \item{Correction}{
  A string containing one of the possible corrections: \code{"None"} (no correction), \code{"ChaoShen"}, \code{"Grassberger"}, \code{"Holste"}, \code{"Bonachela"} or \code{"Best"}, the default value.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
  \item{x}{
  An object to be tested or plotted
  }
  \item{object}{
  A \code{PhyloDiversity} object to be summarized.
  }
  \item{\dots}{
  Additional arguments to be passed to the generic methods.
  }
}
\details{
  The phylogenetic entropy is its generalization of HCDT entropy to unequal species distances (Pavoine et al., 2009).
  
  Diversity is obtained by transforming generalized entropy.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcPhyloDiversity} and choose the \code{Correction}.

  The functions are designed to be used as simply as possible. \code{PhyloDiversity} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcPhyloDiversity} is called. Explicit calls to \code{bcPhyloDiversity} (with bias correction) or to \code{PhyloDiversity.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\note{
  The tree must contain all species of the probability vector. If it contains extra species, computation time will just be increased.
}
\value{
  An object of class \code{PhyloDiversity} is a list:
  \item{Distribution}{The distribution used to calculate diversity}
  \item{Function}{The function used to calculate diversity}
  \item{Tree}{The functional or phylogenetic tree used to calculate diversity}
  \item{Normalized}{Logical. Indicates whether phylodiversity is normalized or proportional to the height of the tree.}
  \item{Type}{The type of diversity ("alpha", "beta" or "gamma").}
  \item{Order}{The order of diversity \eqn{q}.}
  \item{Cuts}{A named vector containing values of neutral diversity along the tree. Names are cut ends, \emph{i.e.} the ends of intervals (the first interval starts at 0 for leaves, the max value is the height of the tree).}
  \item{Total}{A value equal the total diversity (obtained by transforming the total normalized entropy), multiplied by the tree height if \code{Normalize} is FALSE.}
  
  \code{is.PhyloDiversity} returns \code{TRUE} if the object is of class \code{PhyloDiversity}.
  
  \code{summary.PhyloDiversity} returns a summary of the object's value.

  \code{PhyloDiversity} objects can be plotted by \code{\link{plot.PhyloValue}} because \code{PhyloDiversity} objects are also of class \code{PhyloValue}.
}
\references{
  Chao, A., Chiu, C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B} 365(1558): 3599-609.
  
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.
  
  Pavoine, S., Love, M. S. and Bonsall, M. B. (2009). Hierarchical partitioning of evolutionary and ecological patterns in the organization of phylogenetically-structured species assemblages: Application to rockfish (genus: Sebastes) in the Southern California Bight. \emph{Ecology Letters} 12(9): 898-908.
}
\seealso{
  \code{\link{PhyloEntropy}}, \code{\link{Diversity}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ps)
# Calculate the phylogenetic Shannon diversity of the plot
summary(PhyloDiversity(Ps, 1, Paracou618.Taxonomy) -> d)
plot(d)

# Ns is the vector of abundances of the metacommunity
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Calculate the phylogenetic Shannon diversity of the plot
summary(bcPhyloDiversity(Ns, 1, Paracou618.Taxonomy, Correction = "Best") -> d)
plot(d)
}
