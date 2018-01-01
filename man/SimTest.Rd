\name{SimTest}
\alias{SimTest}
\alias{as.SimTest}
\alias{is.SimTest}
\alias{plot.SimTest}
\alias{summary.SimTest}
\title{
  SimTest class
}
\description{
  Methods for objects of type "\code{SimTest}", used to test a value against its distribution under a simulated null hypothesis.
}
\usage{
as.SimTest(RealValue, SimulatedValues)
is.SimTest(x)
\method{plot}{SimTest}(x, Quantiles = c(0.025, 0.975), \dots,
  colValue = "red",lwdValue = 2, ltyValue = 2,
  colQuantiles = "black", lwdQuantiles = 1, ltyQuantiles = 2)
\method{summary}{SimTest}(object, Quantiles = c(0.025, 0.975), \dots)
}
\arguments{
  \item{RealValue}{
  A numeric Value (the actual one).
  }
  \item{SimulatedValues}{
  A numeric vector containing the simulated values.
  }
  \item{x}{
  An object to be tested or plotted.
  }
  \item{object}{
  A \code{SimTest} object to be summarized.
  }
  \item{Quantiles}{
  A vector containing the quantiles of interest.
  }
  \item{colValue}{
  The color of the line representing the real value on the plot.
  }
  \item{lwdValue}{
  The width of the line representing the real value on the plot.
  }
  \item{ltyValue}{
  The line type of the line representing the real value on the plot.
  }
  \item{colQuantiles}{
  The color of the lines representing the quantiles on the plot.
  }
  \item{lwdQuantiles}{
  The width of the lines representing the quantiles on the plot.
  }
  \item{ltyQuantiles}{
  The line type of the lines representing the quantiles on the plot.
  }
  \item{\dots}{
  Additional arguments to be passed to the generic methods.   
  }
}
\details{
  Simulated values should be obtained by simulation. The actual value is compared to simulated quantiles. \code{SimTest} objects can be plotted and summarized.
}
\value{
  \code{SimTest} objects are lists containing:
  \item{RealValue}{The value to test.}
  \item{SimulatedValues}{A vector of simulated values, whose quantiles will be used for the test.}

  \code{is.SimTest} returns \code{TRUE} if the object is of class \code{SimTest}.
  
  \code{summary.SimTest} returns a summary of the object, including the empirical quantile of the real value in the simulated distributon.

}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\examples{
# Set the value to test
Real <- 0.8
# Is it a realization of a Gaussian distribution?
Sims <- rnorm(1000)
# Make a Simtest object
st <- as.SimTest(Real, Sims)
plot(st)
summary(st)
}
