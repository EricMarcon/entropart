\name{rCommunity}
\alias{rCommunity}
\title{
  Random Communities
}
\description{
  Draws random communities according to a probability distribution.
}
\usage{
rCommunity(n, size = sum(NorP), NorP = 1, BootstrapMethod = "Chao2015", S = 300, 
    Distribution = "lnorm",  sd = 1, prob = 0.1, alpha = 40,
    CheckArguments = TRUE)}
\arguments{
  \item{n}{
  The number of communities to draw.
  }
  \item{size}{
  The number of individuals to draw in each community.
  }
  \item{BootstrapMethod}{
  The method used to obtain the probabilities to generate bootstrapped communities from observed abundances. If \code{"Marcon"}, the probabilities are simply the abundances divided by the total number of individuals (Marcon \emph{et al.}, 2012). If \code{"Chao2013"} or \code{"Chao2015"} (by default), a more sophisticated approach is used (see \code{\link{as.ProbaVector}}) following Chao \emph{et al.} (2013) or Chao \emph{et al.} (2015).
  }
  \item{NorP}{
  A numeric vector. Contains either abundances or probabilities.
  }
  \item{S}{
  The number of species.
  }
  \item{Distribution}{
  The distribution of species frequencies. May be \code{"lnorm"} (log-normal), \code{"lseries"} (log-series), \code{"geom"} (geometric) or \code{"bstick"} (broken stick).
  }
  \item{sd}{
  The simulated distribution standard deviation. For the log-normal distribution, this is the standard deviation on the log scale.
  }
  \item{prob}{
  The proportion of ressources taken by successive species.
  }
  \item{alpha}{
  Fisher's alpha.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Communities of fixed \code{size} are drawn in a multinomial distribution according to the distribution of probabilities provided by \code{NorP}.
  
  An abundance vector may be used instead of probabilities, then \code{size} is by default the total number of individuals in the vector. Random communities are built by drawing in a multinomial law following Marcon \emph{et al.} (2012), or trying to estimate the distribution of the actual community with \code{\link{as.ProbaVector}}. If \code{BootstrapMethod =  "Chao2013"}, the distribution is estimated by a single parameter model and unobserved species are given equal probabilities. If \code{BootstrapMethod =  "Chao2015"}, a two-parameter model is used and unobserved species follow a geometric distribution.
  
  Alternatively, the probabilities may be drawn following a classical distribution: either a lognormal (\code{"lnorm"}) one (Preston, 1948) with given standard deviation (\code{sd}; note that the mean is actually a normalizing constant. Its values is set equal to 0 for the simulation of the normal distribution of unnormalized log-abundances), a log-series (\code{"lseries"}) one (Fisher \emph{et al.}, 1943) with parameter \code{alpha}, a geometric (\code{"geom"}) one (Motomura, 1932) with parameter \code{prob}, or a broken stick (\code{"bstick"}) one (MacArthur, 1957). The number of simulated species is fixed by \code{S}, except for \code{"lseries"} where it is obtained from \code{alpha} and \code{size}: \eqn{S=\alpha \ln(1 + \frac{size}{\alpha})}.
  
  Log-normal, log-series and broken-stick distributions are stochastic. The geometric distribution is completely determined by its parameters.
}
\value{
  A vector of species abundances (\code{\link{AbdVector}}) if a single community has been drawn, or a \code{\link{MetaCommunity}} containing simulated communities.
}
\references{
  Chao, A., Wang, Y. T. and Jost, L. (2013). Entropy and the species accumulation curve: a novel entropy estimator via discovery rates of new species. \emph{Methods in Ecology and Evolution} 4(11): 1091-1100.

  Chao, A., Hsieh, T. C., Chazdon, R. L., Colwell, R. K., Gotelli, N. J. (2015) Unveiling the Species-Rank Abundance Distribution by Generalizing Good-Turing Sample Coverage Theory. \emph{Ecology} 96(5): 1189-1201.

  Fisher R.A., Corbet A.S., Williams C.B. (1943) The Relation Between the Number of Species and the Number of Individuals in a Random Sample of an Animal Population. \emph{Journal of Animal Ecology} 12: 42-58.
  
  MacArthur R.H. (1957) On the Relative Abundance of Bird Species. \emph{PNAS} 43(3): 293-295.

  Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The Decomposition of Shannon's Entropy and a Confidence Interval for Beta Diversity. \emph{Oikos} 121(4): 516-522.
  
  Motomura I. (1932) On the statistical treatment of communities. \emph{Zoological Magazine} 44: 379-383.
  
  Preston, F.W. (1948). The commonness, and rarity, of species. \emph{Ecology} 29(3): 254-283.
  
  Reese G. C., Wilson K. R., Flather C. H. (2013) Program SimAssem: Software for simulating species assemblages and estimating species richness. \emph{Methods in Ecology and  Evolution} 4: 891-896.
}
\seealso{
  \code{\link{SpeciesDistribution}} and the program \code{SimAssem} (Reese \emph{et al.}, 2013; not an R package) for more distributions.
}
\examples{
# Generate communities made of 100000 individuals among 300 species and fit them
par(mfrow = c(2,2))
for (d in c("lnorm", "lseries", "geom", "bstick")) {  
rCommunity(n = 1, size = 1E5, S = 300, Distribution = d) -> AbdVec
plot(AbdVec, Distribution = d, main = d)
}
}
