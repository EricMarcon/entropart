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
Simpson(NorP, \dots)
bcSimpson(Ns, Correction = "Best", CheckArguments = TRUE)
\method{Simpson}{ProbaVector}(NorP, \dots, CheckArguments = TRUE, 
  Ps = NULL)
\method{Simpson}{AbdVector}(NorP, Correction="Best", Level = NULL, \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{Simpson}{integer}(NorP, Correction="Best", Level = NULL, \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{Simpson}{numeric}(NorP, Correction="Best", Level = NULL, \dots, 
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
  \item{Correction}{
  A string containing one of the possible corrections accepted by \code{\link{bcTsallis}} or \code{"Lande"}. \code{"Best"}, the default value, is currently \code{"Jackknife"}. Ignored by interpolation and extrapolation.
  }
  \item{Level}{
  The level of interpolation or extrapolation. It may be an a chosen sample size (an integer) or a sample coverage (a number between 0 and 1).
  }
  \item{\dots}{
  Additional arguments. Unused.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Lande's correction has been derived (Lande, 1996; Good, 1953) especially for Simpson entropy, while other corrections are for generalized Tsallis entropy. It is identical to the unbiased estimator proposed by Simpson, although arguments were different. It equals the plug-in etimator multiplied by n/(n-1) where n is the total number of individuals.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}.

  The functions are designed to be used as simply as possible. \code{Simpson} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcSimpson} is called.

  Entropy can be estimated at a specified level of interpolation or extrapolation, either a chosen sample size or sample coverage (Chao et al., 2014), rather than its asymptotic value. Simpson's extrapolated entropy estimator does not rely on the estimation of the asymptotic distribution.
}
\value{
  A named number equal to the calculated entropy. The name is that of the bias correction used.
}
\references{
  Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K., Ellison, A. M (2014). Rarefaction and extrapolation with Hill numbers: A framework for sampling and estimation in species diversity studies.  \emph{Ecological Monographs}, 84(1): 45-67.

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