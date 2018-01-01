\name{TsallisBeta}
\alias{TsallisBeta}
\alias{bcTsallisBeta}
\alias{TsallisBeta.ProbaVector}
\alias{TsallisBeta.AbdVector}
\alias{TsallisBeta.integer}
\alias{TsallisBeta.numeric}
\title{
  Tsallis beta entropy of a community
}
\description{
  Calculates the Tsallis beta entropy of order \eqn{q} of a community belonging to a metacommunity.
}
\usage{
TsallisBeta(NorP, NorPexp = NULL, q = 1, \dots)
bcTsallisBeta(Ns, Nexp = NULL, q, Correction = "Best", CheckArguments = TRUE)
\method{TsallisBeta}{ProbaVector}(NorP, NorPexp = NULL, q = 1, \dots, 
  CheckArguments = TRUE, Ps = NULL, Pexp = NULL)
\method{TsallisBeta}{AbdVector}(NorP, NorPexp = NULL, q = 1, Correction = "Best", \dots, 
  CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
\method{TsallisBeta}{integer}(NorP, NorPexp = NULL, q = 1, Correction = "Best", \dots, 
  CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
\method{TsallisBeta}{numeric}(NorP, NorPexp = NULL, q = 1, Correction = "Best", \dots, 
  CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
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
  A number: the order of entropy. Default is 1 for Shannon entropy.
  }
  \item{Correction}{
  A string containing one of the possible corrections: currently, only \code{"ChaoShen"} or \code{"None"}. \code{"Best"} is the default value, it is equivalent to \code{"ChaoShen"}.
  }
  \item{\dots}{
  Additional arguments. Unused.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The derivation of Tsallis beta entropy can be found in Marcon et al. (2014).
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcTsallisBeta} and choose the \code{Correction}.

  Note that beta entropy value is related to alpha entropy (if \eqn{q} is not 1) and cannot be compared accross communities (Jost, 2007). Beta entropy of a community is not meaningful in general, do rather calculate the \code{\link{BetaDiversity}} of the metacommunity.

  The functions are designed to be used as simply as possible. \code{TsallisBeta} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcTsallisBeta} is called. Explicit calls to \code{bcTsallisBeta} (with bias correction) or to \code{TsallisBeta.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
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
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- Paracou618.MC$Ps
# Probability distribution of the first plot
Ps1 <- Paracou618.MC$Psi[, 1]
# Divergence of order 2 between plot 1 and the whole forest
TsallisBeta(Ps1, Ps, 2)

# Ns is the vector of abundances of the metacommunity
Ns <- Paracou618.MC$Ns
# Abundances in the first plot
Ns1 <- Paracou618.MC$Nsi[, 1]
# Divergence of order 2 between plot 1 and the whole forest, with bias correction
bcTsallisBeta(Ns1, Ns, 2)
}
