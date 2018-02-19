\name{RAC}
\alias{RAClnorm}
\alias{RACgeom}
\alias{RAClseries}
\alias{RACbstick}
\title{
  Fit Distributions to Well-Known Rank Abundance Curves.
}
\description{
  Observed distributions are fitted to classical RAC's.
}
\usage{
RAClnorm(Ns, CheckArguments = TRUE)
RACgeom(Ns, CheckArguments = TRUE)
RAClseries(Ns, CheckArguments = TRUE)
RACbstick(Ns, CheckArguments = TRUE)

}
\arguments{
  \item{Ns}{
  A numeric vector containing species abundances.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  \code{SpeciesDistribution} or integer vectors can be used to fit classical rank-abundance curves (RAC) of classical distributions: \code{"RAClnorm"} for log-normal (Preston, 1948), \code{"RAClseries"} for log-series (Fisher \emph{et al.}, 1943), \code{"RACgeom"} for geometric (Motomura, 1932) or \code{"RACbstick"} for broken stick (MacArthur, 1957). 
  method returns the estimated parameters of the fitted distribution. The broken stick has no parameter, so the maximum abundance is returned.
}
\note{
  Fisher's alpha is estimated to fit the log-series distribution. The estimation is done by the \code{\link{fisher.alpha}} function of package \code{vegan}. It may differ substantially from the estimation returned by \code{optimal.theta} from package \code{untb}.
}
\value{
  A list (the parameters of distributions are returned only if the distribution has been fit):
  \item{Rank}{A numeric vector. The ranks of species in the fitted RAC.}
  \item{Abundance}{The abundance of species in the fitted RAC.}
  \item{mu}{The expectation of the log-normal distribution}
  \item{sigma}{The standard deviation of the log-normal distribution}
  \item{alpha}{Fisher's alpha in the log-series distribution}
  \item{prob}{The proportion of ressources taken by successive species in the geometric distribution}
  \item{max}{The maximum abundance in the broken-stick distribution}
}

\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>, Bruno Herault <Bruno.Herault@ecofog.gf>
}
\references{
  Fisher R.A., Corbet A.S., Williams C.B. (1943) The Relation Between the Number of Species and the Number of Individuals in a Random Sample of an Animal Population. \emph{Journal of Animal Ecology} 12: 42-58.

  MacArthur R.H. (1957) On the Relative Abundance of Bird Species. \emph{PNAS} 43(3): 293-295.

  Motomura I. (1932) On the statistical treatment of communities. \emph{Zoological Magazine} 44: 379-383.
  
  Preston, F.W. (1948). The commonness, and rarity, of species. \emph{Ecology} 29(3): 254-283.
}
\seealso{
  \code{\link{rgeom}}, \code{\link{rlnorm}}, \code{\link{rCommunity}}, \code{\link{plot.SpeciesDistribution}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Fitted parameters
RACln <- RAClnorm(Ns)
RACln$mu
RACln$sigma
RACgeom(Ns)$prob
RAClseries(Ns)$alpha
RACbstick(Ns)$max
}
