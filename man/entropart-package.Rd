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

  In the entropart package, individuals of different "species" are counted in several "communities" which may (or not) be agregated to define a "metacommunity". In the metacommunity, the probability to find a species in the weighted average of probabilities in communities. This is a naming convention, which may correspond to plots in a forest inventory or any data organized the same way.

  Basic functions allow computing diversity of a community. Data is simply a vector of probabilities (summing up to 1) or of abundances (integer values that are numbers of individuals). Calculate entropy with functions such as \code{\link{Tsallis}}, \code{\link{Shannon}}, \code{\link{Simpson}}, \code{\link{Hurlbert}} or \code{\link{GenSimpson}} and explicit diversity (i.e. effective number of species) with \code{\link{Diversity}} and others. By default, the best available estimator of diversity will be used, according to the data.
  
  Communities can be simulated by \code{\link{rCommunity}}, explicitely declared as a species distribution (\code{\link{as.AbdVector}} or \code{\link{as.ProbaVector}}), and plotted.
  
  Phylogenetic entropy and diversity can be calculated if a phylogenetic (or functional), ultrametric tree is provided. See \code{\link{PhyloEntropy}}, \code{\link{Rao}} for examples of entropy and \code{\link{PhyloDiversity}} to calculate phylodiversity, with the state-of-the-art estimation-bias correction. Similarity-based diversity is calculated with \code{\link{Dqz}}, based on a similarity matrix.

  The simplest way to import data is to organize it into two text files. The first file should contain abundance data: the first column named \code{Species} for species names, and a column for each community.

  The second file should contain the community weights in two columns. The first one, named \code{Communities} should contain their names and the second one, named \code{Weights}, their weights.
  
  Files can be read and data imported by code such as:
  
\preformatted{
Abundances <- read.csv(file="Abundances.csv", row.names = 1)
Weights <- read.csv(file="Weights.csv")
MC <- MetaCommunity(Abundances, Weights)
}

  The last line of the code calls the \code{\link{MetaCommunity}} function to create an object that will be used by all metacommunity functions, such as \code{\link{DivPart}} (to partition diversity), \code{\link{DivEst}} (to partition diversity and calculate confidence interval of its estimation) or \code{\link{DivProfile}} (to compute diversity profiles).
  
  A full documentation is available in the vignette. Type: \code{vignette("entropart")}. A quick introuction is in \code{vignette("introduction", "entropart")}.

}
\author{
  Eric Marcon,
  Bruno Herault
}
\references{
  Grabchak, M., Marcon, E., Lang, G., and Zhang, Z. (2017). The Generalized Simpson's Entropy is a Measure of Biodiversity. \emph{Plos One}, 12(3): e0173305.

  Marcon, E. (2015) Practical Estimation of Diversity from Abundance Data. \emph{HAL} 01212435: 1-27.
  
  Marcon, E. and Herault, B. (2015). entropart: An R Package to Measure and Partition Diversity. \emph{Journal of Statistical Software}, 67(8): 1-26.
    
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.
  
  Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The Decomposition of Shannon's Entropy and a Confidence Interval for Beta Diversity. \emph{Oikos} 121(4): 516-522.

  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.

  Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of similarity-based diversity and its bias correction. \emph{HAL} hal-00989454(version 3).
}
\keyword{package}
