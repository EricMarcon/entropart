\name{Simpson}
\alias{Simpson}
\alias{bcSimpson}
\alias{Simpson.ProbaVector}
\alias{Simpson.AbdVector}
\alias{Simpson.integer}
\alias{Simpson.numeric}
\title{
  Simpson entropy of a community
}
\description{
  Calculates the Simpson entropy of a probability vector.
}
\usage{
Simpson(NorP, Correction="Lande", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL)
bcSimpson(Ns, Correction = "Lande", CheckArguments = TRUE)
\method{Simpson}{ProbaVector}(NorP, Correction="Lande", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL)
\method{Simpson}{AbdVector}(NorP, Correction="Lande", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL)
\method{Simpson}{integer}(NorP, Correction="Lande", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL)
\method{Simpson}{numeric}(NorP, Correction="Lande", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL) 
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
    A string containing one of the possible corrections accepted by \code{\link{bcTsallis}} or \code{"Lande"}, the default value (equivalent to \code{"Best"}).
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Lande's correction has been derived (Lande, 1996; Good, 1953) especially for Simpson entropy, while other corrections are for generalized Tsallis entropy. It is identical to the unbiased estimator proposed by Simpson, although arguments were different. It equals the plug-in etimator multiplied by n/(n-1) where n is the total number of individuals.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcSimpson} and choose the \code{Correction}.

  The functions are designed to be used as simply as possible. \code{Simpson} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcSimpson} is called. Explicit calls to \code{bcSimpson} (with bias correction) or to \code{Simpson.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A named number equal to the calculated entropy. The name is that of the bias correction used.
}
\references{
  Good, I. J. (1953). On the Population Frequency of Species and the Estimation of Population Parameters. \emph{Biometrika} 40(3/4): 237-264.

  Lande, R. (1996). Statistics and partitioning of species diversity, and similarity among multiple communities. \emph{Oikos} 76: 5-13.

  Simpson, E. H. (1949). Measurement of diversity. \emph{Nature} 163(4148): 688.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Tsallis}}, \code{\link{bcSimpson}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ns is the total number of trees per species
  Ns <- as.AbdVector(Paracou618.MC$Ns)
  # Whittaker plot
  plot(Ns)

  # Calculate an unbiased estimator of Simpson's index of diversity
  Simpson(Ns)
}