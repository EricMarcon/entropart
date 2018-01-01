\name{Optimal.Similarity}
\alias{Optimal.Similarity}
\title{
  Optimal scale parameter to transform a distance matrix into a similarity matrix
}
\description{
  Calculates the scale parameter \eqn{u} that maximizes the variance of the similarity matrix \eqn{exp(-u*DistanceMatrix)}.
}
\usage{
Optimal.Similarity(Distance, CheckArguments = TRUE)
}
\arguments{
\item{Distance}{
  A distance matrix, \emph{i.e.} a square matrix with zeros on its diagonal or a \code{\link{dist}} object.
}
\item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
}
}
\details{
  The similarity matrix used by \code{\link{Dqz}}) can be optimized following Marcon \emph{et al.} (2014) such that the variance of similarities between pairs of species is maximized. See the example.
}
\value{
  A list:
  \item{u}{The optimal scale \eqn{u}.}
  \item{Matrix}{The optimal similarity matrix \code{Z}.}
}
\references{
  Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of similarity-based diversity and its bias correction. \emph{HAL} hal-00989454(version 3).
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Dqz}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Prepare the similarity matrix. The square root of Paracou618.dist is euclidean.
optimal <- Optimal.Similarity(sqrt(Paracou618.dist))
# Optimal scale
optimal$u
# Calculate diversity of order 2
bcDqz(Paracou618.MC$Ns, 2, optimal$Matrix)
}
