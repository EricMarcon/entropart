\name{DivEst}
\alias{DivEst}
\alias{is.DivEst}
\alias{plot.DivEst}
\alias{autoplot.DivEst}
\alias{summary.DivEst}
\title{
  Diversity Estimation of a metacommunity
}
\description{
  Estimates diversity of a metacommunity.
}
\usage{
DivEst(q = 0, MC, Biased = TRUE, Correction = "Best", Tree = NULL, 
  Normalize = TRUE, Z = NULL, Simulations = 100, 
  ShowProgressBar = TRUE, CheckArguments = TRUE)
is.DivEst(x)
\method{plot}{DivEst}(x, \dots, main = NULL, Which = "All",
  Quantiles = c(0.025, 0.975), colValue = "red", lwdValue = 2, ltyValue = 2,
  colQuantiles = "black", lwdQuantiles = 1, ltyQuantiles = 2)
\method{autoplot}{DivEst}(object, \dots, main = NULL, Which = "All",
  labels = NULL, font.label = list(size=11, face="plain"),
          Quantiles = c(0.025, 0.975), colValue = "red", 
          colQuantiles = "black", ltyQuantiles = 2)
\method{summary}{DivEst}(object, \dots)
}
\arguments{
  \item{q}{
  A number: the order of diversity.
  }
  \item{MC}{
  A \code{\link{MetaCommunity}} object.
  }
  \item{Biased}{
  Logical; if \code{FALSE}, a bias correction is appplied.
  }
  \item{Correction}{
  A string containing one of the possible corrections. The correction must be accepted by \code{\link{DivPart}}. \code{"Best"} is the default value.
  }
  \item{Tree}{
  An object of class \code{\link{hclust}}, "phylo" (see \code{\link[ape]{read.tree}}), \code{\link[ade4]{phylog}} or \code{\link{PPtree}}. The tree must be ultrametric.
  }
  \item{Normalize}{
  If \code{TRUE} (default), diversity is not affected by the height of the tree..\cr
  If \code{FALSE}, diversity is proportional to the height of the tree.
  }
  \item{Z}{
  A relatedness matrix, \emph{i.e.} a square matrix whose terms are all positive, strictly positive on the diagonal. Generally, the matrix is a similarity matrix, \emph{i.e.} the diagonal terms equal 1 and other terms are between 0 and 1.
  }
  \item{Simulations}{
  The number of simulations to build confidence intervals.
  }
  \item{ShowProgressBar}{
  If \code{TRUE} (default), a progress bar is shown.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
  \item{x}{
  An object to be tested or plotted.
  }
  \item{main}{
  The title of the plot.
  }
  \item{Which}{
  May be \code{"Alpha"}, \code{"Beta"} or \code{"Gamma"} to respectively plot the metacommunity's alpha, beta or gamma diversity. If \code{"All"} (default), all three plots are shown.
  }
  \item{labels}{
  Vector of labels to be added to multiple plots. \code{"auto"} is the same as \code{c("a", "b", "c", "d)"}.
  }
  \item{font.label}{
  A list of arguments to customize labels. See \code{\link[ggpubr]{ggarrange}}.
  }
  \item{object}{
  A \code{MCdiversity} object to be summarized or plotted.
  }
  \item{Quantiles}{
  A vector containing the quantiles of interest.
  }
  \item{colValue}{
  The color of the line representing the real value on the plot.
  }
  \item{lwdValue}{
  The width of the line representing the real value on the plot.
  }
  \item{ltyValue}{
  The line type of the line representing the real value on the plot.
  }
  \item{colQuantiles}{
  The color of the lines representing the quantiles on the plot.
  }
  \item{lwdQuantiles}{
  The width of the lines representing the quantiles on the plot.
  }
  \item{ltyQuantiles}{
  The line type of the lines representing the quantiles on the plot.
  }
  \item{\dots}{
  Additional arguments to be passed to the generic methods.
  }
}
\details{
  \code{Divest} estimates the diversity of the metacommunity and partitions it into alpha and beta components.
  
  If \code{Tree} is provided, the phylogenetic diversity is calculated else if \code{Z} is not \code{NULL}, then similarity-based entropy is calculated.
  
  Bootstrap confidence intervals are calculated by drawing simulated communities from a multinomial distribution following the observed frequencies (Marcon et al, 2012; 2014).
}
\value{
  A \code{Divest} object which is a \code{\link{DivPart}} object with an additional item in its list:
  \item{SimulatedDiversity}{A matrix containing the simulated values of alpha, beta and gamma diversity.}

  \code{Divest} objects can be summarized and plotted.
}
\references{
  Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The Decomposition of Shannon's Entropy and a Confidence Interval for Beta Diversity. \emph{Oikos} 121(4): 516-522.
  
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.

  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.
}
\author{
  Eric Marcon <Eric.Marcon@agroparistech.fr>, Bruno Herault <Bruno.Herault@cirad.fr>
}
\seealso{
  \code{\link{DivPart}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Estimate Shannon diversity.
Estimation <- DivEst(q = 1, Paracou618.MC, Biased = FALSE, Correction = "UnveilJ", 
Simulations = 20)
plot(Estimation)
summary(Estimation)
}
