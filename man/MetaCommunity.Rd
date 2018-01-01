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

  The second file should contain the community weights in two columns. The first one, named \code{Communities} should contain their names and the second one, named \code{Weights}, their weights.
  
  Files can be read and data imported by code such as:
  
\preformatted{
Abundances <- read.csv(file="Abundances.csv", row.names = 1)
Weights <- read.csv(file="Weights.csv")
MC <- MetaCommunity(Abundances, Weights)
}

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
