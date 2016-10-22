\name{entropart-package}
\alias{entropart-package}
\alias{entropart}
\docType{package}
\title{
  Entropy Partitioning to Measure Diversity
}
\description{
  Functions to calculate alpha, beta and gamma diversity of communities, including phylogenetic and functional diversity.
  
  Estimation-bias corrections are available.
}
\details{

  In the entropart package, individuals of different "species" are counted in several "communities" which are agregated to define a "metacommunity". In the metacommunity, the probability to find a species in the weighted average of probabilities in communities.
  
  This is a naming convention, which may correspond to plots in a forest inventory or any data organized the same way.
  
  

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
  Abundances <- read.csv(file="Abundances.csv")
  Weights <- read.csv(file="Weights.csv")
  MC <- MetaCommunity(Abundances, Weights)
  }
  
  The last line of the code calls the \code{\link{MetaCommunity}} function to create an object that will be used by all metacommunity functions, 
  such as \code{\link{DivPart}} (to partition diversity), \code{\link{DivEst}} (to partition diversity and calculate confidence interval of its estimation) 
  or \code{\link{DivProfile}} (to compute diversity profiles).
  
  The package also provides community functions to calculate entropy or diversity with bias correction from a vector of probabilities or abundances, such as \code{\link{bcShannon}} (to calculate the reduced-bias estimator of Shannon entropy of a community) or \code{\link{Diversity}} (to calculate the true diversity of a community).
  
  Phylogenetic entropy and diversity can be calculated if a phylogenetic (or functional), ultrametric tree is provided. See \code{\link{PhyloEntropy}}, \code{\link{Rao}} for examples of entropy and \code{\link{bcPhyloDiversity}} to calculate phylodiversity, with bias correction. Similarity-based diversity is calculated with \code{\link{Dqz}} or (with bias correction) \code{\link{bcDqz}}, based on a similarity matrix.

}
\author{
  Eric Marcon,
  Bruno Herault
}
\references{
  Marcon, E. and Herault, B. (2015). entropart: An R Package to Measure and Partition Diversity. \emph{Journal of Statistical Software}, 67(8): 1-26.
    
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.
  
  Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The Decomposition of Shannon's Entropy and a Confidence Interval for Beta Diversity. \emph{Oikos} 121(4): 516-522.

  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.

  Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of similarity-based diversity and its bias correction. \emph{HAL} hal-00989454(version 1).
}
\keyword{ package }
