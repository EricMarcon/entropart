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
  
  Distances are not represented the same way in different formats of trees so choices must be made.
  From version 1.4.8 up, all dendrograms in entropart follow conventions of the package \code{\link{ape}}.
  Edge lengths are those of the \code{\link{phylo}} object, stored in its component \code{$edge.length}. 
  When converted to \code{\link{hclust}}, edge lengths are multiplied by 2 to obtain the values of \code{$heights}, see \code{\link{as.hclust.phylo}}.
  In \code{\link{phylog}} objects, \code{$Wdist} contains values that are twice the square root of distances between tips.
}

\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
