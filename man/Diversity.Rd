\name{Diversity}
\alias{Diversity}
\alias{bcDiversity}
\alias{Diversity.ProbaVector}
\alias{Diversity.AbdVector}
\alias{Diversity.integer}
\alias{Diversity.numeric}
\title{
  HCDT diversity of a community
}
\description{
  Calculates the HCDT (generalized) diversity of order \eqn{q} of a probability vector.
}
\usage{
Diversity(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL)
bcDiversity(Ns, q = 1, Correction = "Best", CheckArguments = TRUE)
\method{Diversity}{ProbaVector}(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL)
\method{Diversity}{AbdVector}(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL)
\method{Diversity}{integer}(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, 
  Ps = NULL, Ns = NULL)
\method{Diversity}{numeric}(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, 
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
  \item{q}{
    A number: the order of diversity. Default is 1.
  }
  \item{Correction}{
    A string containing one of the possible corrections: \code{"None"} (no correction), \code{"ChaoShen"}, \code{"GenCov"}, \code{"Grassberger"}, \code{"Holste"}, \code{"Bonachela"}, \code{"ZhangGrabchak"}, or \code{"ChaoWangJost"}, \code{"Marcon"}, \code{"UnveilC"}, \code{"UnveiliC"}, \code{"UnveilJ"} or \code{"Best"}, the default value.  Currently, \code{"Best"} is \code{"ChaoWangJost"}.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  \code{Diversity} calls \code{\link{Tsallis}} to calculate entropy and transforms it into diversity by calculating its deformed exponential.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcDiversity} and choose the \code{Correction}.
  
  The functions are designed to be used as simply as possible. \code{Diversity} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcDiversity} is called. Explicit calls to \code{bcDiversity} (with bias correction) or to \code{Diversity.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A named number equal to the calculated diversity. The name is that of the bias correction used.
}
\references{
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
  Diversity(Ps, 1)
  # Calculate it with estimation bias correction
  Diversity(Ns, 1)
}
