\name{Richness}
\alias{bcRichness}
\alias{Richness}
\alias{Richness.ProbaVector}
\alias{Richness.AbdVector}
\alias{Richness.integer}
\alias{Richness.numeric}
\title{
  Number of species of a community
}
\description{
  Calculates the number of species from probability vector.  The name is that of the estimator (the bias correction) used.
}
\usage{
Richness(NorP, \dots)
bcRichness(Ns, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, 
              CheckArguments = TRUE)
\method{Richness}{ProbaVector}(NorP, \dots, CheckArguments = TRUE, Ps = NULL)
\method{Richness}{AbdVector}(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, 
  \dots, CheckArguments = TRUE, Ns = NULL)
\method{Richness}{integer}(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, 
  \dots, CheckArguments = TRUE, Ns = NULL)
\method{Richness}{numeric}(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, 
  \dots, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  \item{Correction}{
  A string containing one of the possible corrections: \code{"None"} (no correction), \code{"Jackknife"}, \code{"iChao1"}, or \code{"Chao1"}, the default value.
  }
  \item{Alpha}{
  The risk level, 5\% by default, used to optimize the jackknife order.
  }
  \item{JackOver}{
  If \code{TRUE}, retain the jackknife order immediately superior to the optimal one, usually resulting in the overestimation of the number of species. Default is \code{FALSE}.
  }
  \item{\dots}{
  Additional arguments. Unused.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Bias correction requires the number of individuals. Use \code{bcRichness} and choose the \code{Correction}.
  
  Chao correction techniques are from Chao (1984) and Chiu \emph{et al.} (2015). The Jackknife estimator is calculating by a straight adaptation of the code by Ji-Ping Wang (\code{\link[SPECIES]{jackknife}} in package \code{SPECIES}). The optimal order is selected according to Burnham and Overton (1978; 1979). The argument \code{JackOver} allows selecting one order over the optimal one.
  Many other estimators are available elsewhere, the ones implemented here are necessary for other entropy estimations.

  The functions are designed to be used as simply as possible. \code{Richness} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcRichness} is called. Explicit calls to \code{bcRichness} (with bias correction) or to \code{Richness.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A named number equal to the estimated number of species.
}
\references{
  Burnham, K. P., and Overton,W. S. (1978), Estimation of the Size of a Closed Population When Capture Probabilities Vary Among Animals. \emph{Biometrika}, 65:625-633.

  Burnham, K. P., and Overton,W. S. (1979), Robust Estimation of Population Size When Capture Probabilities Vary Among Animals. \emph{Ecology} 60:927-936.

  Chao, A. (1984) Nonparametric estimation of the number of classes in a population. \emph{Scandinavian Journal of Statistics} 11:265-270.
  
  Chiu, C.-H., Wang, Y.-T., Walther, B. A., Chao, A. (2014) An Improved Nonparametric Lower Bound of Species Richness via a Modified Good-Turing Frequency Formula. \emph{Biometrics} 70(3):671-682.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
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

# Number of observed species
Richness(Ps)
# Estimate the actual number of species
bcRichness(Ns, Correction = "Chao1")
bcRichness(Ns, Correction = "iChao1")
bcRichness(Ns, Correction = "Jackknife")
bcRichness(Ns, Correction = "Jackknife", JackOver=TRUE)
}