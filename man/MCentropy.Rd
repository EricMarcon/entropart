\name{MCentropy}
\alias{MCentropy}
\alias{is.MCentropy}
\alias{plot.MCentropy}
\alias{autoplot.MCentropy}
\alias{summary.MCentropy}
\title{
  Meta-Community entropy class.
}
\description{
  Methods for objects of type "\code{MCentropy}".
}
\usage{
is.MCentropy(x)
\method{plot}{MCentropy}(x, \dots)
\method{autoplot}{MCentropy}(object, \dots)
\method{summary}{MCentropy}(object, \dots)
}
\arguments{
  \item{x}{
  An object to be tested or plotted.
  }
  \item{object}{
  A \code{MCentropy} object to be summarized or plotted.
  }
  \item{\dots}{
  Additional arguments to be passed to the generic methods.
  }
}
\value{
  Meta-community entropy objects are lists containing:
  \item{MetaCommunity}{The name of the \code{\link{MetaCommunity}} object containing inventory data.}
  \item{Method}{The method used to calculate entropy ("HCDT", "Similarity-based").}
  \item{Type}{The type of entropy ("alpha", "beta" or "gamma").}
  \item{Order}{The order of entropy \eqn{q}.}
  \item{Correction}{The estimation bias correction used to calculate entropy.}
  \item{Tree}{The phylogenetic or functional tree used to calculate phyloentropy.}
  \item{Normalized}{Logical. Indicates whether phyloentropy is normalized or proportional to the height of the tree.}
  \item{Z}{The matrix used to calculate similarity-based entropy.}
  \item{Weights}{A vector containing the weights of communities.}
  \item{Communities}{A vector containing the entropy of communities.}
  \item{Total}{The total entropy.}
  
  \code{is.MCentropy} returns \code{TRUE} if the object is of class \code{MCentropy}.

  \code{summary.MCentropy} returns a summary of the object's value.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
