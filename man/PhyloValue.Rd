\name{PhyloValue}
\alias{PhyloValue}
\alias{is.PhyloValue}
\alias{plot.PhyloValue}
\alias{summary.PhyloValue}
\title{
  Phylogenetic entropy of diversity.
}
\description{
  Entropy or diversity against the height of the phylogenetic or functional tree.
}
\usage{
is.PhyloValue(x)
\method{plot}{PhyloValue}(x, xlab = expression(italic("T")), ylab = NULL, main = NULL, \dots)
\method{summary}{PhyloValue}(object, \dots)
}
\arguments{
  \item{x}{
    An object of class \code{PhyloValue}, including \code{\link{PhyloDiversity}} and \code{\link{PhyloEntropy}} objects.
  }
  \item{xlab}{
    The X axis label, "T" by default for Time.
  }
  \item{ylab}{
    The Y axis label. if \code{NULL} (by default), "Entropy" or "Diversity" or nothing is chosen according to the object class.
  }
  \item{main}{
    The main title of the plot. if \code{NULL} (by default), a default value is used.
  }
  \item{object}{
    A \code{PhyloValue} object to be summarized.
  }
  \item{\dots}{
    Additional arguments to be passed to \code{\link{plot}}.
  }
}
\details{
  \code{PhyloValue} objects are the result of \code{\link{PhyloApply}}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
