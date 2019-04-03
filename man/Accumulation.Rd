\name{Accumulation}
\alias{DivAC}
\alias{EntAC}
\alias{as.AccumCurve}
\alias{is.AccumCurve}
\alias{autoplot.AccumCurve}
\alias{plot.AccumCurve}

\title{
  Diversity accumulation.
}
\description{
  Diversity and Entropy Accumulation Curves represent the accumulation of entropy with respect to the sample size.
}
\usage{
as.AccumCurve(x, y, low = NULL, high = NULL)
is.AccumCurve(x)
EntAC(Ns, q = 0, n.seq = 1:sum(Ns), PCorrection="Chao2015", Unveiling="geom",
  RCorrection="Rarefy", NumberOfSimulations = 0, Alpha = 0.05, 
  ShowProgressBar = TRUE, CheckArguments = TRUE)
DivAC(Ns, q = 0, n.seq = 1:sum(Ns), PCorrection="Chao2015", Unveiling="geom",
  RCorrection="Rarefy", NumberOfSimulations = 0, Alpha = 0.05, 
  ShowProgressBar = TRUE, CheckArguments = TRUE)
\method{plot}{AccumCurve}(x, ..., main = NULL, 
         xlab = "Sample Size", ylab = NULL, ylim = NULL,
         LineWidth = 2, ShadeColor = "grey75", BorderColor = "red")
\method{autoplot}{AccumCurve}(object, ..., main = NULL, 
         xlab = "Sample Size", ylab = NULL, 
         ShadeColor = "grey75", alpha = 0.3, BorderColor = "red")}
\arguments{
  \item{x}{
  An object. A numeric vector in \code{as.AccumCurve}.
  }
  \item{object}{
  An object.
  }
  \item{y}{
  A numeric vector.
  }
  \item{low}{
  A numeric vector.
  }
  \item{high}{
  A numeric vector.
  }
  \item{Ns}{
  A numeric vector containing species abundances.
  }
  \item{q}{
  A number: the order of diversity. Default is 1.
  }
  \item{n.seq}{
  A sequence of numbers. Accumulation will be calculated at each value.
  }
  \item{PCorrection}{
  A string containing one of the possible corrections to estimate a probability distribution in \code{\link{as.ProbaVector}}: \code{"Chao2015"} is the default value. Used only for extrapolation and \code{q} different from 0, 1, 2.
  }
  \item{Unveiling}{
  A string containing one of the possible unveiling methods to estimate the probabilities of the unobserved species in \code{\link{as.ProbaVector}}: \code{"geom"} (geometric: the unobserved species distribution is geometric) is the default value. Used only for extrapolation and \code{q} different from 0, 1, 2.
  }
  \item{RCorrection}{
  A string containing a correction recognized by \code{\link{Richness}} to evaluate the total number of species in \code{\link{as.ProbaVector}}. \code{"Rarefy"} is the default value to estimate the number of species such that the entropy of the asymptotic distribution rarefied to the observed sample size equals the observed entropy of the data. Used only for extrapolation and \code{q} different from 0, 1, 2. If \code{q} is 0 (extrapolation of richness), \code{"Rarefy"} is taken for \code{"Jackknife"}.
  }
  \item{NumberOfSimulations}{
  The number of Simulations to build confidence intervals.
  }
  \item{Alpha}{
  The risk level, 5\% by default.
  }
  \item{\dots}{
  Additional arguments to be passed to \code{\link{plot}}. Unused elsewhere.
  }
  \item{main}{
  The main title of the plot. if \code{NULL} (by default), there is no title.
  }
  \item{xlab}{
  The X axis label, "Rank" by default.
  }
  \item{ylab}{
  The Y axis label. if \code{NULL} (by default), "Probability" or "Abundance" is chosen according to the object class.
  }
  \item{ylim}{
  The interval of y values plotted.
  }
  \item{LineWidth}{
  The width of the line that represents the actual profile.
  }
  \item{ShadeColor}{
  The color of the shaded confidence envelope.
  }
  \item{BorderColor}{
  The color of the bounds of the confidence envelope.
  }
  \item{alpha}{
  Opacity of the confidence enveloppe, between 0 and 1.
  }
  \item{ShowProgressBar}{
  If \code{TRUE} (default), a progress bar is shown.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  \code{DivAC} or \code{EntAC} estimate the diversity or entropy accumulation curve of a distribution.
  See \code{\link{Tsallis}} for details about the computation of entropy at each level of interpolation and extrapolation.
  In accumulation curves, extrapolation if done by estimating the asymptotic distribution of the community and estimating entropy at different levels by interpolation. The asymptotic richess is adjusted so that the extrapolated part of the accumulation joins the observed value at the sample size.

  \code{AccumCurve} objects include \code{EntAC} and \code{DivAC} objects for entropy and diversity accumulation.
They generalize the classical Species Accumulation Curves (SAC) which are diversity accumulation of order $q=0$.
  
  \code{as.AccumCurve} transforms two vectors (where \code{x} is the sammple size and \code{y} the accumulation) into an object of class \code{AccumCurve}.

  \code{AccumCurve} objects can be plotted with either \code{plot} or \code{autoplot methods}.
}
\value{
  A \code{DivAC} or an \code{EntAC} object. Both are \code{AccumCurve} objects, which are a list:
  \item{x}{The sample size.}
  \item{y}{The value of entropy or diversity.}
  \item{low}{The lower bound of the confidence envelope of the estimation.}
  \item{high}{The upper bound of the confidence envelope of the estimation.}
 
  Attibutes \code{"Size"} and \code{"Value"} contain the actual sample size and the corresponding diversity or entropy.
  
  \code{AccumCurve} objects can be summarized and plotted.
}
\references{
  Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K., Ellison, A. M (2014). Rarefaction and extrapolation with Hill numbers: A framework for sampling and estimation in species diversity studies.  \emph{Ecological Monographs}, 84(1): 45-67.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>, Bruno Herault <Bruno.Herault@ecofog.gf>
}
\seealso{
  \code{\link{Tsallis}}, \code{\link{Diversity}}
}

\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Accumulation curve of Simpson's diversity
autoplot(DivAC(Ns, q=2))
}
