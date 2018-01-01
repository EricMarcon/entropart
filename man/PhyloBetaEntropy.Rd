\name{PhyloBetaEntropy}
\alias{PhyloBetaEntropy}
\alias{bcPhyloBetaEntropy}
\alias{PhyloBetaEntropy.ProbaVector}
\alias{PhyloBetaEntropy.AbdVector}
\alias{PhyloBetaEntropy.integer}
\alias{PhyloBetaEntropy.numeric}
\title{
  Phylogenetic Beta Entropy of a community
}
\description{
  Calculates the phylogenetic beta entropy of order \eqn{q} of a a community belonging to a metacommunity.
}
\usage{
PhyloBetaEntropy(NorP, NorPexp = NULL, q = 1, Tree, Normalize = TRUE, \dots)
bcPhyloBetaEntropy(Ns, Nexp, q = 1, Tree, Normalize = TRUE, 
  Correction = "Best", CheckArguments = TRUE)
\method{PhyloBetaEntropy}{ProbaVector}(NorP, NorPexp = NULL, q = 1, Tree, Normalize = TRUE, 
  \dots, CheckArguments = TRUE, Ps = NULL, Pexp = NULL)
\method{PhyloBetaEntropy}{AbdVector}(NorP, NorPexp = NULL, q = 1, Tree, Normalize = TRUE, 
  Correction = "Best", \dots, CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
\method{PhyloBetaEntropy}{integer}(NorP, NorPexp = NULL, q = 1, Tree, Normalize = TRUE, 
  Correction = "Best", \dots, CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
\method{PhyloBetaEntropy}{numeric}(NorP, NorPexp = NULL, q = 1, Tree, Normalize = TRUE, 
  Correction = "Best", \dots, CheckArguments = TRUE, Ps = NULL, Ns = NULL, 
  Pexp = NULL, Nexp = NULL) 
}
\arguments{
  \item{Ps}{
  The probability vector of species of the community.
  }
  \item{Pexp}{
  The probability vector of species of the metacommunity.
  }
  \item{Ns}{
  A numeric vector containing species abundances of the community.
  }
  \item{Nexp}{
  A numeric vector containing species abundances of the metacommunity.
  }
  \item{NorP}{
  A numeric vector, an integer vector, an abundance vector (\code{\link{AbdVector}}) or a probability vector (\code{\link{ProbaVector}}). Contains either abundances or probabilities of the community.
  }
  \item{NorPexp}{
  A numeric vector, an integer vector, an abundance vector (\code{\link{AbdVector}}) or a probability vector (\code{\link{ProbaVector}}). Contains either abundances or probabilities of the metacommunity.
  }
  \item{q}{
  A number: the order of entropy. Default is 1.
  }
  \item{Tree}{
  An object of class \code{\link{hclust}}, \code{\link{phylo}}, \code{\link{phylog}} or \code{\link{PPtree}}. The tree must be ultrametric.
  }
  \item{Normalize}{
  If \code{TRUE} (default), the entropy returned by the function is normalized by the height of the tree (it is the weighted average value of the entropy in each slice).\cr
  If \code{FALSE}, it is the unnormalized weighted sum of the results.
  }
  \item{Correction}{
  A string containing one of the possible corrections: currently, only \code{"ChaoShen"}. \code{"Best"} is the default value, it is equivalent to \code{"ChaoShen"}.
  }
  \item{\dots}{
  Additional arguments. Unused.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The phylogenetic entropy is the generalization of HCDT entropy to unequal species distances (Pavoine et al., 2009).
  
  Calculation relies on \code{\link{TsallisBeta}} and \code{\link{PhyloApply}}.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcPhyloBetaEntropy} and choose the \code{Correction}.
  
  Note that beta entropy value is related to alpha entropy (if \eqn{q} is not 1) and cannot be compared accross communities (Jost, 2007). Beta entropy of a community is not meaningful in general, do rather calculate the \code{\link{PhyloDiversity}} of the metacommunity.

  The functions are designed to be used as simply as possible. \code{PhyloBetaEntropy} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcPhyloBetaEntropy} is called. Explicit calls to \code{bcPhyloBetaEntropy} (with bias correction) or to \code{PhyloBetaEntropy.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A \code{\link{PhyloEntropy}} object containing entropy values at each cut of the tree.
}
\references{
  Jost (2007), Partitioning diversity into independent alpha and beta components. \emph{Ecology} 88(10): 2427-2439.

  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.
  
  Pavoine, S., Love, M. S. and Bonsall, M. B. (2009). Hierarchical partitioning of evolutionary and ecological patterns in the organization of phylogenetically-structured species assemblages: Application to rockfish (genus: Sebastes) in the Southern California Bight. \emph{Ecology Letters} 12(9): 898-908.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{TsallisBeta}}, \code{\link{bcPhyloBetaEntropy}}, \code{\link{PhyloDiversity}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ps)
# Probability distribution of the first plot
Ps1 <- as.ProbaVector(Paracou618.MC$Psi[, 1])
# Calculate the phylogenetic Shannon beta entropy of the plot
summary(PhyloBetaEntropy(Ps1, Ps, 1, Paracou618.Taxonomy) -> e)
plot(e)

# Ns is the vector of abundances of the metacommunity
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Abundances in the first plot
Ns1 <- as.AbdVector(Paracou618.MC$Nsi[, 1])
# Calculate the phylogenetic Shannon beta entropy of the plot
summary(bcPhyloBetaEntropy(Ns1, Ns, 1, Paracou618.Taxonomy, Correction = "Best") -> e)
plot(e)
}
