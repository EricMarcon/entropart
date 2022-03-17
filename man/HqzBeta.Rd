\name{HqzBeta}
\alias{HqzBeta}
\alias{bcHqzBeta}
\alias{HqzBeta.ProbaVector}
\alias{HqzBeta.AbdVector}
\alias{HqzBeta.integer}
\alias{HqzBeta.numeric}
\title{
  Similarity-based beta entropy of a community
}
\description{
  Calculates the similarity-based beta entropy of order \eqn{q} of a community belonging to a metacommunity.
}
\usage{
HqzBeta(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), \dots)
bcHqzBeta(Ns, Nexp = NULL, q = 1, Z = diag(length(Ns)), Correction = "Best",
          CheckArguments = TRUE)
\method{HqzBeta}{ProbaVector}(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)),
  \dots, CheckArguments = TRUE, Ps = NULL, Pexp = NULL)
\method{HqzBeta}{AbdVector}(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), Correction = "Best", 
  \dots, CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
\method{HqzBeta}{integer}(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), Correction = "Best", 
  \dots, CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
\method{HqzBeta}{numeric}(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), Correction = "Best", 
  \dots, CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
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
  A number, the order of diversity. Default is 1.
  }
  \item{Z}{
  A relatedness matrix, \emph{i.e.} a square matrix whose terms are all positive, strictly positive on the diagonal. Generally, the matrix is a similarity matrix, \emph{i.e.} the diagonal terms equal 1 and other terms are between 0 and 1. Default is the identity matrix to calculate neutral entropy.
  }
  \item{Correction}{
  A string containing one of the possible corrections: currently, no correction is available so \code{"Best"},  the default value, is equivalent to \code{"None"}.
  }
  \item{\dots}{
  Additional arguments. Unused.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The derivation of similarity-based beta entropy can be found in Marcon \emph{et al.} (2014).
  
  Bias correction requires the number of individuals.

  Note that beta entropy value is related to alpha entropy (if \eqn{q} is not 1) and cannot be compared accross communities (Jost, 2007). Beta entropy of a community is not meaningful in general, do rather calculate the \code{\link{BetaDiversity}} of the metacommunity.
  
  The functions are designed to be used as simply as possible. \code{HqzBeta} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcHqzBeta} is called. Explicit calls to \code{bcHqzBeta} (with bias correction) or to \code{HqzBeta.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A named number equal to the calculated entropy. The name is that of the bias correction used.
}
\references{
  Jost (2007), Partitioning diversity into independent alpha and beta components. \emph{Ecology} 88(10): 2427-2439.

  Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of similarity-based diversity and its bias correction. \emph{HAL} hal-00989454(version 3).
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ps)
# Probability distribution of the first plot
Ps1 <- as.ProbaVector(Paracou618.MC$Psi[, 1])
# Prepare the similarity matrix
DistanceMatrix <- as.matrix(Paracou618.dist)
# Similarity can be 1 minus normalized distances between species
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
# Divergence of order 2 between plot 1 and the whole forest
HqzBeta(Ps1, Ps, q=2, Z)
}
