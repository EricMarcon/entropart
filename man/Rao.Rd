\name{Rao}
\alias{Rao}
\alias{bcRao}
\alias{Rao.ProbaVector}
\alias{Rao.AbdVector}
\alias{Rao.integer}
\alias{Rao.numeric}
\title{
  Rao Quadratic Entropy of a Community
}
\description{
  Calculates Rao's quadratic entropy of a community described by a probability vector and a phylogenetic / functional tree.
}
\usage{
Rao(NorP, Tree, \dots)
bcRao(Ns, Tree, Correction="Lande", CheckArguments = TRUE)
\method{Rao}{ProbaVector}(NorP, Tree, \dots, CheckArguments = TRUE, Ps = NULL)
\method{Rao}{AbdVector}(NorP, Tree, Correction = "Lande", \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{Rao}{integer}(NorP, Tree, Correction = "Lande", \dots, 
  CheckArguments = TRUE, Ns = NULL)
\method{Rao}{numeric}(NorP, Tree, Correction = "Lande", \dots, 
  CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
}
\arguments{
  \item{Ps}{
  A probability vector, summing to 1.
  }
  \item{Ns}{
  A numeric vector containing species abundances.
  }
  \item{NorP}{
  A numeric vector, an integer vector, an abundance vector (\code{\link{AbdVector}}) or a probability vector (\code{\link{ProbaVector}}). Contains either abundances or probabilities.
  }
  \item{Tree}{
  An object of class \code{\link{hclust}}, \code{\link{phylo}}, \code{\link[ade4]{phylog}} or \code{\link{PPtree}}. The tree must be ultrametric.
  }
  \item{Correction}{
  A string containing one of the possible corrections accepted by \code{\link{bcTsallis}} or \code{"Lande"}, the default value (equivalent to \code{"Best"}).
  }
  \item{\dots}{
  Additional arguments. Unused.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Bias correction requires the number of individuals. Use \code{bcRao} and choose the \code{Correction}.
    
  The unbiased estimator of Rao's entropy is identical to that of Simpson's entropy because Rao's entropy is a linear sum of Simson entropies, all of them calculated from the same number of individuals (Marcon and Herault, 2014). It equals the plug-in etimator multiplied by n/(n-1) where n is the total number of individuals.

  The functions are designed to be used as simply as possible. \code{Tsallis} is a generic method. If its first argument is an abundance vector, an integer vector or a numeric vector which does not sum to 1, the bias corrected function \code{bcTsallis} is called. Explicit calls to \code{bcTsallis} (with bias correction) or to \code{Tsallis.ProbaVector} (without correction) are possible to avoid ambiguity. The \code{.integer} and \code{.numeric} methods accept \code{Ps} or \code{Ns} arguments instead of \code{NorP} for backward compatibility.
}
\value{
  A named number equal to the calculated entropy. The name is that of the bias correction used.
}
\references{
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.

  Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified approach. \emph{Theoretical Population Biology} 21: 24-43.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{bcPhyloDiversity}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Species probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ns)
# Calculate Rao's quadratic entropy of the plot
Rao(Ps, Paracou618.Taxonomy)
}
