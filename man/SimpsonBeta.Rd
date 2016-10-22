\name{SimpsonBeta}
\alias{SimpsonBeta}
\alias{bcSimpsonBeta}
\alias{SimpsonBeta.ProbaVector}
\alias{SimpsonBeta.AbdVector}
\alias{SimpsonBeta.integer}
\alias{SimpsonBeta.numeric}
\title{
  Simpson beta entropy of a community
}
\description{
  Calculates the Simpson beta entropy of a community belonging to a metacommunity.
}
\usage{
SimpsonBeta(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL)
bcSimpsonBeta(Ns, Nexp, Correction = "Best", CheckArguments = TRUE)
\method{SimpsonBeta}{ProbaVector}(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL)
\method{SimpsonBeta}{AbdVector}(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL)
\method{SimpsonBeta}{integer}(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL)
\method{SimpsonBeta}{numeric}(NorP, NorPexp = NULL, Correction = "Best", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
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
  \item{Correction}{
    A string containing one of the possible corrections: currently, only \code{"ChaoShen"}, identical to \code{"Best"}.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The derivation of Tsallis beta entropy (Simpson is Tsallis of order 2) can be found in Marcon et al. (2014).
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcSimpsonBeta} and choose the \code{Correction}.
  
  Note that Simpson beta entropy value is related to Simpson alpha entropy value and cannot be compared accross communities (Jost, 2007). Beta entropy of a community is not meaningful in general, do rather calculate the \code{\link{BetaDiversity}} of order 2 of the metacommunity.

  The functions are designed to be used as simply as possible. \code{SimpsonBeta} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcSimpsonBeta} is called. Explicit calls to \code{bcSimpsonBeta} (with bias correction) or to \code{SimpsonBeta.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A named number equal to the calculated entropy. The name is that of the bias correction used.
}
\references{
  Jost (2007), Partitioning diversity into independent alpha and beta components. \emph{Ecology} 88(10): 2427-2439.
  
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Simpson}}, \code{\link{bcSimpsonBeta}}, \code{\link{BetaDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- as.ProbaVector(Paracou618.MC$Ps)
  # Probability distribution of the first plot
  Ps1 <- as.ProbaVector(Paracou618.MC$Psi[, 1])

  # Simpson beta entropy of the plot
  SimpsonBeta(Ps1, Ps)
  # Transform into diversity
  expq(SimpsonBeta(Ps1, Ps)/(1-Simpson(Ps1)), 2)
  
  # Ns is the vector of abundances of the metacommunity
  Ns <- as.AbdVector(Paracou618.MC$Ns)
  # Abundances in the first plot
  Ns1 <- as.AbdVector(Paracou618.MC$Nsi[, 1])
  # Reduced-bias Shannon beta entropy of the plot
  bcSimpsonBeta(Ns1, Ns)
}
