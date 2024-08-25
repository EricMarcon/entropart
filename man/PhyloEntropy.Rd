\name{PhyloEntropy}
\alias{PhyloEntropy}
\alias{bcPhyloEntropy}
\alias{PhyloEntropy.ProbaVector}
\alias{PhyloEntropy.AbdVector}
\alias{PhyloEntropy.integer}
\alias{PhyloEntropy.numeric}
\alias{is.PhyloEntropy}
\alias{summary.PhyloEntropy}
\title{
  Phylogenetic Entropy of a community
}
\description{
  Calculates the phylogenetic entropy of order \eqn{q} of a probability vector.
}
\usage{
PhyloEntropy(NorP, q = 1, Tree, Normalize = TRUE, \dots)
bcPhyloEntropy(Ns, q = 1, Tree, Normalize = TRUE, Correction = "Best", 
  SampleCoverage = NULL, CheckArguments = TRUE)
\method{PhyloEntropy}{ProbaVector}(NorP, q = 1, Tree, Normalize = TRUE, 
  \dots, CheckArguments = TRUE, Ps = NULL)
\method{PhyloEntropy}{AbdVector}(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", 
  \dots, CheckArguments = TRUE, Ns = NULL)
\method{PhyloEntropy}{integer}(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", 
  \dots, CheckArguments = TRUE, Ns = NULL)
\method{PhyloEntropy}{numeric}(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", 
  \dots, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
is.PhyloEntropy(x)
\method{summary}{PhyloEntropy}(object, \dots)
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
  A number: the order of entropy. Default is 1.
  }
  \item{Tree}{
  An object of class \code{\link{hclust}}, "phylo" (see \code{\link[ape]{read.tree}}), \code{\link[ade4]{phylog}} or \code{\link{PPtree}}. The tree must be ultrametric.
  }
  \item{Normalize}{
  If \code{TRUE} (default), the \code{Total} entropy returned by the function is normalized by the height of the tree (it is the weighted average value of the entropy in each slice).\cr
  If \code{FALSE}, it is the unnormalized weighted sum of the results.
  }
  \item{Correction}{
  A string containing one of the possible corrections supported by \code{\link{Tsallis}}.
  }
  \item{SampleCoverage}{
  The sample coverage of \code{Ns} calculated elsewhere. Used to calculate the gamma diversity of meta-communities, see details. 
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
  \item{x}{
  An object to be tested or plotted
  }
  \item{object}{
  A \code{PhyloEntropy} object to be summarized.
  }
  \item{\dots}{
  Additional arguments to be passed to the generic methods.
  }
}
\details{
  The phylogenetic entropy is its generalization of HCDT entropy to unequal species distances (Pavoine et al., 2009).
  
  Calculation relies on \code{\link{Tsallis}} and \code{\link{PhyloApply}}.
  
  Intervals separate two cuts in a tree: no node is found at heights contained in an interval.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcPhyloEntropy} and choose the \code{Correction}.

  The functions are designed to be used as simply as possible. \code{PhyloEntropy} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcPhyloEntropy} is called. Explicit calls to \code{bcPhyloEntropy} (with bias correction) or to \code{PhyloEntropy.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.

  The size of a metacommunity (see \code{\link{MetaCommunity}}) is unknown so it has to be set according to a rule which does not ensure that its abundances are integer values. Then, classical bias-correction methods do not apply. Providing the \code{SampleCoverage} argument allows applying the \code{"ChaoShen"} and \code{"Grassberger"} corrections to estimate quite well the entropy. \code{\link{DivPart}} and \code{\link{GammaEntropy}} functions use this tweak.
}
\note{
  The tree must contain all species of the probability vector. If it contains extra species, computation time will just be increased.
}
\value{
  An object of class \code{PhyloEntropy} is a list:
  \item{Distribution}{The distribution used to calculate entropy}
  \item{Function}{The function used to calculate entropy}
  \item{Tree}{The functional or phylogenetic tree used to calculate entropy}
  \item{Normalized}{Logical. Indicates whether phyloentropy is normalized or proportional to the height of the tree.}
  \item{Type}{The type of entropy ("alpha", "beta" or "gamma").}
  \item{Order}{The order of entropy \eqn{q}.}
  \item{Cuts}{A named vector containing values of neutral entropy along the tree. Names are cut ends, \emph{i.e.} the ends of intervals (the first interval starts at 0 for leaves, the max value is the height of the tree).}
  \item{Total}{A value equal the total entropy multiplied by the tree height if \code{Normalize} is FALSE.}
  
  \code{is.PhyloEntropy} returns \code{TRUE} if the object is of class \code{PhyloEntropy}.
  
  \code{summary.PhyloEntropy} returns a summary of the object's value.

  \code{PhyloEntropy} objects can be plotted by \code{\link{plot.PhyloValue}} because \code{PhyloEntropy} objects are also of class \code{PhyloValue}.
}
\references{
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.
  
  Pavoine, S., Love, M. S. and Bonsall, M. B. (2009). Hierarchical partitioning of evolutionary and ecological patterns in the organization of phylogenetically-structured species assemblages: Application to rockfish (genus: Sebastes) in the Southern California Bight. \emph{Ecology Letters} 12(9): 898-908.
}
\seealso{
  \code{\link{Tsallis}}, \code{\link{PhyloDiversity}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ps)
# Calculate the phylogenetic Shannon entropy of the plot
summary(PhyloEntropy(Ps, 1, Paracou618.Taxonomy) -> e)
plot(e)

# Ns is the vector of abundances of the metacommunity
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Calculate the phylogenetic Shannon entropy of the plot
summary(bcPhyloEntropy(Ns, 1, Paracou618.Taxonomy, Correction = "Best") -> e)
plot(e)
}
