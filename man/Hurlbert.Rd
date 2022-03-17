\name{Hurlbert}
\alias{Hurlbert}
\alias{bcHurlbert}
\alias{Hurlbert.ProbaVector}
\alias{Hurlbert.AbdVector}
\alias{Hurlbert.integer}
\alias{Hurlbert.numeric}
\alias{HurlbertD}
\alias{bcHurlbertD}
\alias{HurlbertD.ProbaVector}
\alias{HurlbertD.AbdVector}
\alias{HurlbertD.integer}
\alias{HurlbertD.numeric}
\title{
  Hurlbert's Index and Explicit Diversity
}
\description{
  Calculates the Hurlbert entropy of order \eqn{k} of a probability or abundance vector, and its effective number of species.
}
\usage{
Hurlbert(NorP, k = 2, \dots)
bcHurlbert(Ns, k = 2, CheckArguments = TRUE)
\method{Hurlbert}{ProbaVector}(NorP, k = 2, \dots, 
  CheckArguments = TRUE, Ps = NULL)
\method{Hurlbert}{AbdVector}(NorP, k = 2, \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{Hurlbert}{integer}(NorP, k = 2, \dots,
  CheckArguments = TRUE, Ns = NULL)
\method{Hurlbert}{numeric}(NorP, k = 2, \dots,
  CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
HurlbertD(NorP, k = 2, ...)
bcHurlbertD(Ns, k = 2, CheckArguments = TRUE)
\method{HurlbertD}{ProbaVector}(NorP, k = 2, \dots,
  CheckArguments = TRUE, Ps = NULL)
\method{HurlbertD}{AbdVector}(NorP, k = 2, \dots,
  CheckArguments = TRUE, Ns = NULL)
\method{HurlbertD}{integer}(NorP, k = 2, \dots,
  CheckArguments = TRUE, Ns = NULL)
\method{HurlbertD}{numeric}(NorP, k = 2, \dots,
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
  \item{k}{
  A number: the order of diversity. Default is 2 for Simpson's diversity.
  }
  \item{\dots}{
  Additional arguments. Unused.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Hurlbert's index of diversity (1971) of order \eqn{k} is the expected number of species in a sample of size \eqn{k}.
  
  Bias correction requires the number of individuals. Use \code{bcHurlbert}. It is limited to orders \eqn{k} less than or equal to the number of individuals in the community.

  The effective number of species \code{HurlbertD} (explicit diversity) has been derived by Dauby & Hardy (2012). It is calculated numerically. \code{bcHurlbertD} calculates it from the bias-corrected index \code{bcHurlbert}.

  The functions are designed to be used as simply as possible. \code{Hurlbert} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcHurlbert} is called. Explicit calls to \code{bcHurlbert} (with bias correction) or to \code{Hurlbert.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A named number equal to the calculated index or diversity. The name is either "Biased" or "Unbiased", depending on the estimator used.
}
\references{
  Dauby G. & Hardy O.J. (2012) Sampled-based estimation of diversity sensu stricto by transforming Hurlbert diversities into effective number of species. \emph{Ecography} 35(7): 661-672.

  Hurlbert (1971) The Nonconcept of Species Diversity: A Critique and Alternative Parameters. \emph{Ecology} 52(4): 577-586.
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

# Calculate Hurlbert entropy of order 2, equal to Simpson's index of diversity
Hurlbert(Ps, 2) 
# Calculate an unbiased estimator of Hurlbert entropy of order 2
Hurlbert(Ns, 2)
}
