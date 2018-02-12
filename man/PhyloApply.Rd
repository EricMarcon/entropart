\name{PhyloApply}
\alias{PhyloApply}
\title{
  Apply a Function over a Phylogenetic Tree
}
\description{
  Cuts the tree into slices separated by nodes, applies the function to each slice and returns the weighted (by slice lengths) sum of the results.
}
\usage{
PhyloApply(Tree, FUN, NorP, Normalize = TRUE, dfArgs = NULL,
           \dots, CheckArguments = TRUE)
}
\arguments{
  \item{Tree}{
  An object of class \code{\link{hclust}}, \code{\link{phylo}}, \code{\link{phylog}} or \code{\link{PPtree}}. The tree must be ultrametric.
  }
  \item{FUN}{
  The function to be applied to each interval of the tree. 
  }
  \item{NorP}{
  A numeric vector or a two-column matrix. Contains either abundances or probabilities. Two-column matrices should contain the observed abundances (or probabilities) in the first column and the expected ones in the second column, to allow using beta diversity functions.
  }
  \item{Normalize}{
  If \code{TRUE} (default), the \code{Total} value returned by \code{Function} is normalized by the height of the tree (it is the weighted average value of the result in each slice).\cr
  If \code{FALSE}, it is the unnormalized weighted sum of the results.
  }
  \item{dfArgs}{
  A dataframe. Columns are arguments for \code{FUN}: their names are those of valid arguments. Values will be passed to \code{FUN} in each slice of the tree, starting from the tips. The number of lines must equal the number of slices.
  }
  \item{\dots}{
  Further arguments to pass to \code{Function}.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  This funtion is generally not used directly. It is a tool to calculate \code{\link{PhyloEntropy}} and \code{\link{PhyloDiversity}}.
  
  Intervals (slices) separate two cuts (nodes) in a tree: no node is found at heights contained in an interval.
  
  Objects of class \strong{PPtree} are returned by \code{\link{Preprocess.Tree}}.
  
  \dots allow passing arguments to the function but they can't change along the tree. If necessary, \code{dfArgs} allow passing a different value for each slice of the tree.
}
\value{
  An object of class \code{PhyloValue}. It is a list:
  \item{Distribution}{The distribution used to calculate the value}
  \item{Function}{The function used to calculate the value}
  \item{Tree}{The functional or phylogenetic tree used to calculate the value}
  \item{Normalized}{Logical. Indicates whether phylovalue is normalized or proportional to the height of the tree.}
  \item{Cuts}{A named vector containing values along the tree. Names are cut ends, \emph{i.e.} the ends of intervals (the first interval starts at 0 for leaves, the max value is the height of the tree).}
  \item{Corrections}{A named vector containing the correction used by FUN to obtain each value of \code{Cuts}. Names are those of \code{Cuts}.}
  \item{Total}{The total value, multiplied by the tree height if \code{Normalize} is FALSE.}
}
\references{
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Preprocess.Tree}}
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Plot the taxonomy
plot(Paracou618.Taxonomy, type="fan", show.tip.label=FALSE)
# Calculate the mean number of trees (individuals) per species 
# (Cuts are 1=species, 2=genus, 3=family)
summary(PhyloApply(Paracou618.Taxonomy, mean, Paracou618.MC$Ns, TRUE))
}
