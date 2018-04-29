\name{AbdFreqCount}
\alias{AbdFreqCount}
\title{
  Abundance Frequency Count of a Community
}
\description{
  Counts the number of species observed the same number of times.
}
\usage{
AbdFreqCount(Ns, Level = NULL, Estimator = "Best", CheckArguments = TRUE)
}
\arguments{
  \item{Ns}{
  A numeric vector containing species abundances.
  }
  \item{Level}{
  The level of interpolation or extrapolation. It may be an a chosen sample size (an integer) or a sample coverage (a number between 0 and 1).
  }
  \item{Estimator}{
  A string containing one of the possible estimators (ignored if \code{Level} is \code{NULL}: \code{"UnveilJ"} only at the moment. \code{"Best"} is for \code{"UnveilJ"}.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The Abundance Frequency Count (Chao \emph{et al.}, 2015) is the number of species observed each number of times.
  
  It is a way to summarize the species distribution.
  
  It can be estimated at a specified level of interpolation or extrapolation. Extrapolation relies on the estimation of eq. (5) of Chao et al. (2014). The \code{"UnveilJ"} estimator builds the full distribution (Chao et al., 2015) of probabilities.
}
\value{
  A two-column matrix. The first column contains the number of observations, the second one the number of species observed this number of times.
}
\references{
  Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K., and  Ellison, A. M. (2014). Rarefaction and extrapolation with Hill numbers: A framework for sampling and estimation in species diversity studies.  \emph{Ecological Monographs}, 84(1): 45–67.

  Chao, A., Hsieh, T. C., Chazdon, R. L., Colwell, R. K., Gotelli, N. J. (2015) Unveiling the Species-Rank Abundance Distribution by Generalizing Good-Turing Sample Coverage Theory. \emph{Ecology} 96(5): 1189-1201.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{PhyloEntropy}}, \code{\link{ChaoPD}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Ns is the vector of abundances of the first plot
Ns <- Paracou618.MC$Nsi[, 1]

# Return the abundance frequency count
(AbdFreqCount(Ns) -> afc)
plot(afc, xlab="Number of observations", ylab="Number of species")
lines(afc)
}
