\name{GenSimpson}
\alias{GenSimpson}
\alias{bcGenSimpson}
\alias{GenSimpson.ProbaVector}
\alias{GenSimpson.AbdVector}
\alias{GenSimpson.integer}
\alias{GenSimpson.numeric}
\alias{GenSimpsonD}
\alias{bcGenSimpsonD}
\alias{GenSimpsonD.ProbaVector}
\alias{GenSimpsonD.AbdVector}
\alias{GenSimpsonD.integer}
\alias{GenSimpsonD.numeric}
\title{
  Generalized Simpson's Entropy and Diversity
}
\description{
  Calculates the Generalized Simpson's entropy of order \eqn{r} of a probability or abundance vector, and its effective number of species.
}
\usage{
GenSimpson(NorP, r = 1, \dots)
bcGenSimpson(Ns, r = 1, CheckArguments = TRUE)
\method{GenSimpson}{ProbaVector}(NorP, r = 1, \dots, 
  CheckArguments = TRUE, Ps = NULL)
\method{GenSimpson}{AbdVector}(NorP, r = 1, \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{GenSimpson}{integer}(NorP, r = 1, \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{GenSimpson}{numeric}(NorP, r = 1, \dots, 
  CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
GenSimpsonD(NorP, r = 1, \dots)
bcGenSimpsonD(Ns, r = 1, CheckArguments = TRUE)
\method{GenSimpsonD}{ProbaVector}(NorP, r = 1, \dots,
  CheckArguments = TRUE, Ps = NULL)
\method{GenSimpsonD}{AbdVector}(NorP, r = 1, \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{GenSimpsonD}{integer}(NorP, r = 1, \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{GenSimpsonD}{numeric}(NorP, r = 1, \dots, 
  CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  \item{r}{
  A number: the order of diversity. Default is 1 for Simpson's diversity.
  }
  \item{\dots}{
  Additional arguments. Unused.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The Generalized Simpson's Entropy (Zhang and Zhou, 2010) of order \eqn{r} is, in the species accumulation curve, the probability for the individual sampled in rank \eqn{r+1} to belong to a new species. It is a measure of diversity so long as \eqn{r} is lower than the number of species (Grabchak \emph{et al.}, 2016). 
  
  Bias correction requires the number of individuals. Use \code{bcGenSimpson}. It is limited to orders \eqn{r} less than or equal to the number of individuals in the community.

  The effective number of species \code{GenSimpsonD} (explicit diversity) has been derived by Grabchak \emph{et al.} (2016).

  The functions are designed to be used as simply as possible. \code{GenSimpson} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcGenSimpson} is called. Explicit calls to \code{bcGenSimpson} (with bias correction) or to \code{GenSimpson.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A named number equal to the calculated index or diversity. The name is either "Biased" or "Unbiased", depending on the estimator used.
}
\note{
  The unbiased estimator is calculated by the \code{\link{GenSimp.z}} function of the \code{EntropyEstimation} package. 
}
\references{
  Grabchak, M., Marcon, E., Lang, G., and Zhang, Z. (2017). The Generalized Simpson's Entropy is a Measure of Biodiversity. \emph{Plos One}, 12(3): e0173305.

  Zhang Z. and Zhou J. (2010). Re-parameterization of multinomial distributions and diversity indices. \emph{Journal of Statistical Planning and Inference} 140(7): 1731-1738.
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Species probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ns)
# Whittaker plot
plot(Ns)

# Calculate GenSimpson entropy of order 1, equal to Simpson's index of diversity
GenSimpson(Ps, 1) 
# Calculate an unbiased estimator of GenSimpson diversity of order 100
GenSimpsonD(Ns, 100)
}
