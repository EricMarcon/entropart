\name{KLq}
\alias{KLq}
\title{
  Generalized Kullback-Leibler divergence
}
\description{
  Calculates the generalized Kullback-Leibler divergence between an observed and an expected probability distribution.
}
\usage{
KLq(Ps, Pexp, q = 1, CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
  The observed probability vector.
  }
  \item{Pexp}{
  The expected probability vector.
  }
  \item{q}{
  A number: the order of entropy. Default is 1.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The generalized Kullback-Leibler divergence (Borland et al., 1998) converges to the Kullback-Leibler divergence (Kullback and Leibler, 1951) when \eqn{q} tends to 1.
  It is used to calculate the generalized beta entropy (Marcon et al., 2014).
}
\value{
  A number equal to the generalized Kullback-Leibler divergence between the probability distributions.
}
\references{
  Borland, L., Plastino, A. R. and Tsallis, C. (1998). Information gain within nonextensive thermostatistics. \emph{Journal of Mathematical Physics} 39(12): 6490-6501.
  
  Kullback, S. and Leibler, R. A. (1951). On Information and Sufficiency. The Annals of Mathematical \emph{Statistics} 22(1): 79-86.
  
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{TsallisBeta}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- Paracou618.MC$Ps
# Probability distribution of the first plot
Ps1 <- Paracou618.MC$Psi[, 1]
# Divergence of order 2 between the first plot and the whole forest
KLq(Ps1, Ps, 2)
}
