\name{Dqz}
\alias{Dqz}
\alias{bcDqz}
\alias{Dqz.ProbaVector}
\alias{Dqz.AbdVector}
\alias{Dqz.integer}
\alias{Dqz.numeric}
\title{
  Similarity-based diversity of a community
}
\description{
  Calculates the diversity of order \eqn{q} of a probability vector according to a similarity matrix.
}
\usage{
Dqz(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", 
  CheckArguments = TRUE, Ps = NULL, Ns = NULL)
bcDqz(Ns, q = 1, Z = diag(length(Ns)), Correction = "Best", CheckArguments = TRUE)
\method{Dqz}{ProbaVector}(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", 
  CheckArguments = TRUE, Ps = NULL, Ns = NULL)
\method{Dqz}{AbdVector}(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", 
  CheckArguments = TRUE, Ps = NULL, Ns = NULL)
\method{Dqz}{integer}(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", 
  CheckArguments = TRUE, Ps = NULL, Ns = NULL)
\method{Dqz}{numeric}(NorP, q = 1, Z = diag(length(Ps)), Correction = "Best", 
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
  \item{Z}{
    A relatedness matrix, \emph{i.e.} a square matrix whose terms are all positive, strictly positive on the diagonal. Generally, the matrix is a similarity matrix, \emph{i.e.} the diagonal terms equal 1 and other terms are between 0 and 1. Default is the identity matrix to calculate neutral diversity.
  }
  \item{Correction}{
    A string containing one of the possible corrections: \code{"None"} (no correction), \code{"HorvitzThomson"}, \code{"MarconZhang"} or \code{"Best"}, the default value. The \code{"MarconZhang"} correction assumes a similarity matrix.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Diversity is calculated following Leinster and Cobbold (2012): it is the reciprocal of the (generalized) average (of order \code{q}) of the community species ordinariness.
  
  A similarity matrix is used (as for \code{\link{Dqz}}), not a distance matrix as in Ricotta and Szeidl (2006). See the example.

  Bias correction requires the number of individuals. Use \code{bcHqz} and choose the \code{Correction}.
  Correction techniques are from Marcon \emph{et al.} (2014).
  
  Currently, the \code{"Best"} correction is the max value of \code{"HorvitzThomson"} and \code{"MarconZhang"}.

  The functions are designed to be used as simply as possible. \code{Dqz} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcDqz} is called. Explicit calls to \code{bcDqz} (with bias correction) or to \code{Dqz.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A named number equal to the calculated diversity. The name is that of the bias correction used.
}
\references{
  Leinster, T. and Cobbold, C. (2012). Measuring diversity: the importance of species similarity. \emph{Ecology} 93(3): 477-489.

  Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of similarity-based diversity and its bias correction. \emph{HAL} hal-00989454(version 1).
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Hqz}}, \code{\link{PhyloDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Prepare the similarity matrix
  DistanceMatrix <- as.matrix(Paracou618.dist)
  # Similarity can be 1 minus normalized distances between species
  Z <- 1 - DistanceMatrix/max(DistanceMatrix)
  # Calculate diversity of order 2
  Dqz(Paracou618.MC$Ns, 2, Z)
}
