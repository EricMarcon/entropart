\name{Preprocess.MC}
\alias{Preprocess.MC}
\title{
  Preprocessing of metacommunity data
}
\description{
  Calculates summary statistics of a metacommunity
}
\usage{
Preprocess.MC(Nsi, Wi)
}
\arguments{
  \item{Nsi}{
  A matrix containing abundance data, species in line, communities in column.
  }
  \item{Wi}{
  Normalized community weights (sum to 1).
  }
}
\details{
  This is an internal function. It is used by \code{\link{MetaCommunity}} for real data and \code{\link{DivEst}} for simulations.
}
\value{
An object of class \strong{MetaCommunity}. It is a list:
  \item{Nsi}{A matrix containing abundance data, species in line, communities in column}
  \item{Nspecies}{The number of species}
  \item{Ncommunities}{The number of communities}
  \item{Wi}{Normalized community weights (sum to 1)}
  \item{N}{The total number of individuals}
  \item{Ni}{A vector containing the number of individuals of each community}
  \item{Ns}{A vector containing the number of individuals of each species}
  \item{Psi}{A matrix whose columns are the probability vectors of communities (each of them sums to 1)}
  \item{Ps}{A vector containing the probability vector of the metacommunity}
  \item{SampleCoverage}{The sample coverage of the metacommunity}
  \item{SampleCoverage.communities}{A vector containing the sample coverages of each community}
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\keyword{internal}
