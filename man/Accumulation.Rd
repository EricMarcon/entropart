\name{Accumulation}
\alias{EntAC}
\alias{as.AccumCurve}
\alias{is.AccumCurve}
\alias{autoplot.AccumCurve}
\alias{plot.AccumCurve}

\title{
  Diversity accumulation.
}
\description{
  Diversity and entropy Accumulation Curves represent the accumulation of entropy with respect to the sample size.
}
\usage{
as.AccumCurve(x, y, low = NULL, high = NULL)
is.AccumCurve(x)
EntAC(Ns, q = 0, n.seq = 1:sum(Ns), Correction = "Best", 
      Simulations = 0, Alpha = 0.05, CheckArguments = TRUE)
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
  \item{Ns}{
  A numeric vector containing species abundances.
  }
  \item{q}{
  A number: the order of diversity. Default is 1.
  }
  \item{n.seq}{
  A sequence of numbers. Accumulation will be calculated at each value.
  }
  \item{Correction}{
  A string containing one of the possible asymptotic corrections of \code{\link{Tsallis}}: \code{"None"} (no correction) and  \code{"Best"} are two of them.
  }
  \item{Simulations}{
  The number of simulations to build confidence intervals.
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
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  \code{EntAC} estimates the entropy accumulation curve of a distribution.

  \code{AccumCurve} objects include \code{EntAC} and \code{DivAC} objects for entropy and diversity accumulation.
They generalize the classical Species Accumulation Curves (SAC) which are diversity accumulationof order $q=0$.
  
  \code{as.AccumCurve} transforms two vectors (where \code{x} is the sammple size and \code{y} the accumulation) into an object of class \code{AccumCurve}.

  \code{AccumCurve} objects can be plotted with either \code{plot} or \code{autoplot methods}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>, Bruno Herault <Bruno.Herault@ecofog.gf>
}

\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Accumulation curve of Shannon's entropy from 1 to 100 indivuduals
autoplot(EntAC(Ns), n.seq = 1:100)
}
