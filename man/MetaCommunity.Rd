\name{MetaCommunity}
\alias{MetaCommunity}
\alias{is.MetaCommunity}
\alias{plot.MetaCommunity}
\alias{summary.MetaCommunity}
\title{
  Metacommunity class
}
\description{
  Methods for objects of type "\code{MetaCommunity}".
}
\usage{
MetaCommunity(Abundances, Weights = rep(1, ncol(Abundances)))
is.MetaCommunity(x)
\method{summary}{MetaCommunity}(object, \dots)
\method{plot}{MetaCommunity}(x, \dots)
}
\arguments{
  \item{Abundances}{
  A dataframe containing the number of observations (lines are species, columns are communities). The first column of the dataframe may contain the species names.
  }
  \item{Weights}{
  A vector of positive numbers equal to community weights or a dataframe containing a vector named \code{Weights}. It does not have to be normalized. Weights are equal by default.
  }
  \item{x}{
    An object to be tested or plotted.
  }
  \item{object}{
    A \code{MetaCommunity} object to be summarized.
  }
  \item{\dots}{
    Additional arguments to be passed to the generic methods.
  }
}
\details{
  In the entropart package, individuals of different "species" are counted in several "communities" which are agregated to define a "metacommunity".
  
  This is a naming convention, which may correspond to plots in a forest inventory or any data organized the same way.
  
  Alpha and beta entropies of communities are summed according to \code{Weights} and the probability to find a species in the metacommunity is the weighted average of probabilities in communities.
  
  
  The simplest way to import data is to organize it into two text files. The first file should contain abundance data: the first column named \code{Species} for species names, and a column for each community.

  \tabular{lcc}{
    \code{Species}        \tab \code{NameOfCommunity1} \tab \code{NameOfCommunity2}\cr
    \code{NameOfSpecies1} \tab \code{1}                \tab \code{5}               \cr
    \code{NameOfSpecies2} \tab \code{4}                \tab \code{2}               \cr
    \code{...}            \tab \code{...}              \tab \code{...}
  }

  The second file should contain the community weights and be organized as follows:
  
  \tabular{lc}{
    \code{Communities}      \tab \code{Weights}\cr
    \code{NameOfCommunity1} \tab \code{3}      \cr
    \code{NameOfCommunity2} \tab \code{1}
  }

  Files can be read and data imported by code such as:
  
  \preformatted{
  Abundances <- read.csv(file="Abundances.csv", row.names = 1)
  Weights <- read.csv(file="Weights.csv")
  MC <- MetaCommunity(Abundances, Weights)
  }

  The number of individuals of the meta-community \code{$N} is not explicitly specified. It is calculated as the number of individuals of the lower-weighted community divided by its weight. For example, if the community with the smallest weight (5\%) contains 100 individuals, the total number of individuals is set to 2000. The number of individuals per species of the meta-community \code{$Ns} is the product of \code{$N} and \code{$Ps}. In the special case where weights are proportional to the number of individuals of communities, the resulting number of individuals of the meta-community is simply the sum of the number of individuals of communities. In other cases, the numbers of individuals of the meta-community may not be integers. This does not appear to be a problem.
}
\value{
An object of class \strong{MetaCommunity} is a list:
  \item{Nsi}{A matrix containing abundance data, species in line, communities in column.}
  \item{Ns}{A vector containing the number of individuals of each species.}
  \item{Ni}{A vector containing the number of individuals of each community.}
  \item{N}{The total number of individuals.}
  \item{Psi}{A matrix whose columns are the probability vectors of communities (each of them sums to 1).}
  \item{Wi}{A vector containing the normalized community weights (sum to 1).}
  \item{Ps}{A vector containing the probability vector of the metacommunity.}
  \item{Nspecies}{The number of species.}
  \item{Ncommunities}{The number of communities.}
  \item{SampleCoverage}{The sample coverage of the metacommunity.}
  \item{SampleCoverage.communities}{A vector containing the sample coverages of each community.}
  
  \code{is.MetaCommunity} returns \code{TRUE} if the object is of class \code{MetaCommunity}.

  \code{summary.MetaCommunity} returns a summary of the object's value.
  
  \code{plot.MetaCommunity} plots it.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\examples{
  # Use BCI data from vegan package
  if (require(vegan, quietly = TRUE)) {
    # Load BCI data (number of trees per species in each 1-ha plot of a tropical forest)
    data(BCI)
    # BCI dataframe must be transposed (its lines are plots, not species)
    BCI.df <- as.data.frame(t(BCI))
    # Create a metacommunity object from a matrix of abundances and a vector of weights
    # (here, all plots have a weight equal to 1)
    MC <- MetaCommunity(BCI.df)
  }
}
