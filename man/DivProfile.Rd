\name{DivProfile}
\alias{DivProfile}
\alias{is.DivProfile}
\alias{plot.DivProfile}
\alias{summary.DivProfile}
\title{
  Diversity Profile of a metacommunity
}
\description{
  Calculate the diversity profiles (alpha, beta, gamma) of a metacommunity.
}
\usage{
DivProfile(q.seq = seq(0, 2, 0.1), MC, Biased = TRUE, Correction = "Best", 
  Tree = NULL, Normalize = TRUE, Z = NULL, 
  NumberOfSimulations = 0, Alpha = 0.05, CheckArguments = TRUE)
is.DivProfile(x)
\method{plot}{DivProfile}(x, \dots, main = NULL, xlab = "Order of Diversity",
  ylab = NULL, Which = "All", 
  LineWidth = 2, ShadeColor = "grey75", BorderColor = "red")
\method{summary}{DivProfile}(object, \dots)
}
\arguments{
  \item{q.seq}{
  A numeric vector.
  }
  \item{MC}{
  A \code{\link{MetaCommunity}} object.
  }
  \item{Biased}{
  Logical; if \code{FALSE}, a bias correction is appplied.
  }
  \item{Correction}{
  A string containing one of the possible corrections.
	
  The correction must be accepted by \code{\link{AlphaEntropy}}, \code{\link{BetaEntropy}} and \code{\link{GammaEntropy}}. \code{"Best"} is the default value.
  }
  \item{Tree}{
  An object of class \code{\link{hclust}}, \code{\link{phylo}}, \code{\link{phylog}} or \code{\link{PPtree}}. The tree must be ultrametric.
  }
  \item{Normalize}{
  If \code{TRUE} (default), diversity is not affected by the height of the tree.\cr
  If \code{FALSE}, diversity is proportional to the height of the tree.
  }
  \item{Z}{
  A relatedness matrix, \emph{i.e.} a square matrix whose terms are all positive, strictly positive on the diagonal. Generally, the matrix is a similarity matrix, \emph{i.e.} the diagonal terms equal 1 and other terms are between 0 and 1.
  }
  \item{NumberOfSimulations}{
  The number of simulations to run, 0 by default.
  }
  \item{Alpha}{
  The risk level, 5\% by default.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
  \item{x}{
  An object to be tested or plotted.
  }
  \item{main}{
  The main title of the plot. Ignored if \code{Which = "All"}.
  }
  \item{xlab}{
  The x axis label of the plots.
  }
  \item{ylab}{
  The y axis label of the plot. Ignored if \code{Which = "All"}.
  }
  \item{Which}{
  May be \code{"Communities"}, \code{"Alpha"}, \code{"Beta"} or \code{"Gamma"} to respectively plot the alpha diversity of communities or the metacommunity's alpha, beta or gamma diversity. If \code{"All"} (default), all four plots are shown.
  }
  \item{LineWidth}{
  The width of the line that represents the actual profile.
  }
  \item{ShadeColor}{
  The color of the shaded confidence envelope.
  }
  \item{BorderColor}{
  The color of the bounds of the confidence envelope.
  }
  \item{object}{
  A \code{MCdiversity} object to be summarized.
  }
  \item{\dots}{
  Additional arguments to be passed to the generic methods.
  }
}
\details{
  If \code{Tree} is provided, the phylogenetic diversity is calculated.
  
  \code{DivPart} partitions the diversity of the metacommunity into alpha and beta components. It supports estimation-bias correction.
  
  If \code{Tree} is provided, the phylogenetic diversity is calculated else if \code{Z} is not \code{NULL}, then similarity-based entropy is calculated.
  
  Beta diversity/entropy is calculated from Gamma and Alpha when bias correction is required, so community values are not available.
  
  If \code{NumberOfSimulations} is greater than 0, a bootstrap confidence interval is produced by simulating communities from a multinomial distribution following the observed frequencies (Marcon et al, 2012; 2014) and calculating their profiles.

}
\value{
  A \code{DivProfile} object. It is a list:
  \item{MetaCommunity}{The name of the \code{\link{MetaCommunity}} object containing inventory data.}
  \item{Order}{A vector containing the values of \code{q}.}
  \item{Biased}{Logical. If \code{FALSE}, bias corrected values of diversity have been computed.}
  \item{Correction}{The estimation bias correction used to calculate diversity. Usually a string, but it may be a list if different corrections have been used in the estimation of phylodiversity.}
  \item{Method}{The method used to calculate entropy ("HCDT", "Similarity-based").}
  \item{Tree}{The phylogenetic or functional tree used to calculate phylodiversity.}
  \item{Normalized}{Logical. Indicates whether phylodiversity is normalized or proportional to the height of the tree.}
  \item{Z}{The matrix used to calculate similarity-based entropy.}
  \item{CommunityAlphaDiversities}{A matrix containing the alpha diversity of each community.}
  \item{TotalAlphaDiversity}{A vector containing the alpha diversity of communities for each order.}
  \item{BetaDiversity}{A vector containing the beta diversity of communities for each order.}
  \item{GammaDiversity}{A vector containing the gamma diversity of the metacommunity for each order.}
  \item{CommunityAlphaEntropies}{A matrix containing the alpha entropy of each community.}
  \item{TotalAlphaEntropy}{A vector containing the alpha entropy of communities for each order.}
  \item{BetaEntropy}{A vector containing the beta entropy of communities for each order.}
  \item{GammaEntropy}{A vector containing the gamma entropy of the metacommunity for each order.}
  \item{Confidence envelopes}{Total Alpha, Beta and Gamma Entropy and Diversity may come with a confidence envelope whose value is stored in twelve more vectors named suffixed Low or High, such as \code{GammaEntropyLow}}
  
  \code{DivProfile} objects can be summarized and plotted.
}
\references{
  Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The Decomposition of Shannon's Entropy and a Confidence Interval for Beta Diversity. \emph{Oikos} 121(4): 516-522.
  
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.
  
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf> and Bruno Herault <Bruno.Herault@ecofog.gf>
}
\seealso{
  \code{\link{DivPart}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Estimate diversity.
Profile <- DivProfile(q.seq = seq(0, 2, 0.1), Paracou618.MC, Biased = FALSE)
plot(Profile)
summary(Profile)
}
