\name{Diversity}
\alias{Diversity}
\alias{bcDiversity}
\alias{Diversity.ProbaVector}
\alias{Diversity.AbdVector}
\alias{Diversity.integer}
\alias{Diversity.numeric}
\title{
  Hill number of a community
}
\description{
  Calculates the HCDT (generalized) diversity of order \eqn{q} of a probability vector.
}
\usage{
Diversity(NorP, q = 1, \dots)
bcDiversity(Ns, q = 1, Correction = "Best", CheckArguments = TRUE)
\method{Diversity}{ProbaVector}(NorP, q = 1, \dots, 
  CheckArguments = TRUE, Ps = NULL)
\method{Diversity}{AbdVector}(NorP, q = 1, Correction = "Best", Level = NULL,
  PCorrection="Chao2015", Unveiling="geom", RCorrection="Rarefy", \dots,
  CheckArguments = TRUE, Ns = NULL)
\method{Diversity}{integer}(NorP, q = 1, Correction = "Best", Level = NULL,
  PCorrection="Chao2015", Unveiling="geom", RCorrection="Rarefy", \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{Diversity}{numeric}(NorP, q = 1, Correction = "Best", Level = NULL,
  PCorrection="Chao2015", Unveiling="geom", RCorrection="Rarefy", \dots, 
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
  \item{q}{
  A number: the order of diversity. Default is 1.
  }
  \item{Correction}{
  A string containing one of the possible asymptotic estimators: \code{"None"} (no correction), \code{"ChaoShen"}, \code{"GenCov"}, \code{"Grassberger"}, \code{"Holste"}, \code{"Bonachela"}, \code{"ZhangGrabchak"}, or \code{"ChaoJost"}, \code{"Marcon"}, \code{"UnveilC"}, \code{"UnveiliC"}, \code{"UnveilJ"} or \code{"Best"}, the default value. Currently, \code{"Best"} is \code{"UnveilJ"}.
  }
  \item{Level}{
  The level of interpolation or extrapolation. It may be an a chosen sample size (an integer) or a sample coverage (a number between 0 and 1).
  }
  \item{PCorrection}{
  A string containing one of the possible corrections to estimate a probability distribution in \code{\link{as.ProbaVector}}: \code{"Chao2015"} is the default value. Used only for extrapolation.
  }
  \item{Unveiling}{
  A string containing one of the possible unveiling methods to estimate the probabilities of the unobserved species in \code{\link{as.ProbaVector}}: \code{"geom"} (the unobserved species distribution is geometric) is the default value. Used only for extrapolation.
  }
  \item{RCorrection}{
  A string containing a correction recognized by \code{\link{Richness}} to evaluate the total number of species in \code{\link{as.ProbaVector}}. \code{"Rarefy"} is the default value to estimate the number of species such that the diversity of the asymptotic distribution rarefied to the observed sample size equals the observed diversity of the data. Used only for extrapolation.
  }
  \item{\dots}{
  Additional arguments. Unused.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  \code{Diversity} calls \code{\link{Tsallis}} to calculate entropy and transforms it into diversity by calculating its deformed exponential.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. 
  See \code{\link{Tsallis}} for details.
  
  The functions are designed to be used as simply as possible. 
  \code{Diversity} is a generic method. 
  If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcDiversity} is called.

  Diversity can be estimated at a specified level of interpolation or extrapolation, either a chosen sample size or sample coverage (Chao et al., 2014), rather than its asymptotic value.
  See \code{\link{Tsallis}} for details.
}
\value{
  A named number equal to the calculated diversity. The name is that of the bias correction used.
}
\references{
  Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K., Ellison, A. M (2014). Rarefaction and extrapolation with Hill numbers: A framework for sampling and estimation in species diversity studies.  \emph{Ecological Monographs}, 84(1): 45-67.

  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Tsallis}}, \code{\link{expq}}, \code{\link{AbdVector}}, \code{\link{ProbaVector}}
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
# Calculate diversity of order 1, i.e. Shannon's diversity
Diversity(Ps, q=1)
# Calculate it with estimation bias correction (asymptotic estimator)
Diversity(Ns, q=1)
# Extrapolate it up to 99.9\% sample coverage (close to the asymptotic estimator)
Diversity(Ns, q=1, Level=0.999)
# Rarefy it to half the sample size
Diversity(Ns, q=1, Level=sum(Ns)/2)
}
