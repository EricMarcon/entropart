\name{PDFD}
\alias{PDFD}
\title{
  Phylogenetic Diversity / Functional Diversity of a Community
}
\description{
  Calculates Faith's PD / Petchey and Gaston' FD of a community described by a probability vector and a phylogenetic / functional tree.
}
\usage{
PDFD(Ps, Tree, CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
  A probability vector, summing to 1.
  }
  \item{Tree}{
  An object of class \code{\link{hclust}}, \code{\link{phylo}}, \code{\link[ade4]{phylog}} or \code{\link{PPtree}}. The tree must be ultrametric.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  PD and FD are defined as the total legth of the branches of the tree.
  
  The probability vector is used to select branches: branches with probability 0 are eliminated.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}.
  
  Use \code{\link{bcPhyloDiversity}}\code{(Ps, 0, Tree)} and choose the \code{Correction}.
}
\value{
  A named number equal to the calculated diversity. The name is that of the bias correction used.
}
\references{
  Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. \emph{Biological Conservation} 61(1): 1-10.
  
  Petchey, O. L. and Gaston, K. J. (2002). Functional diversity (FD), species richness and community composition. \emph{Ecology Letters} 5: 402-411.
}
\seealso{
  \code{\link{bcPhyloDiversity}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- Paracou618.MC$Ps
# Calculate the phylogenetic Shannon diversity of the plot
PDFD(Ps, Paracou618.Taxonomy)
}
