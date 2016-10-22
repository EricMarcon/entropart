\name{PPtree}
\alias{PPtree}
\alias{is.PPtree}
\alias{plot.PPtree}
\title{
  Preprocessed Trees.
}
\description{
  Methods for objects of type "\code{PPtree}".
}
\usage{
is.PPtree(x)
\method{plot}{PPtree}(x, \dots)
}
\arguments{
  \item{x}{
    An object to be tested or plotted
  }
  \item{\dots}{
    Additional arguments to be passed to the generic methods.
  }
}
\value{
  An object of class \strong{PPtree} is a list:
  \item{phyTree}{A \code{\link{phylo}} tree}
  \item{hTree}{A \code{\link{hclust}} tree}
  \item{Height}{The height of the tree, that is to say the distance between root and leaves}
  \item{Cuts}{A vector. Cut times of the tree (the distance from nodes to leaves)}
  \item{Intervals}{A vector. The lengths of intervals between cuts}

  \code{is.PPtree} returns \code{TRUE} if the object is of class \code{PPtree}.
  
  \code{plot.PPtree} plots it.
}
\note{
  Versions up to 1.3 contained a \code{\link{phylog}} tree, now deprecated in \code{ade4}. A \code{\link{phylo}} tree is now used.
}

\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
