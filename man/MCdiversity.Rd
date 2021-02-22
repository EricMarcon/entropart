\name{MCdiversity}
\alias{MCdiversity}
\alias{is.MCdiversity}
\alias{plot.MCdiversity}
\alias{autoplot.MCdiversity}
\alias{summary.MCdiversity}
\title{
  Meta-Community diversity class.
}
\description{
  Methods for objects of type "\code{MCdiversity}".
}
\usage{
is.MCdiversity(x)
\method{plot}{MCdiversity}(x, \dots)
\method{autoplot}{MCdiversity}(object, col = ggplot2::GeomCol$default_aes$fill,
          border = ggplot2::GeomCol$default_aes$colour, \dots)
\method{summary}{MCdiversity}(object, \dots)
}
\arguments{
  \item{x}{
  An object to be tested or plotted.
  }
  \item{object}{
  A \code{MCdiversity} object to be summarized or plotted.
  }
  \item{col}{
  The color used to fill the bars. See "Color Specification" in \code{\link{par}}.
  }
  \item{border}{
  The color of the borders around the bars. See \code{\link{hist}}.
  }
  \item{\dots}{
  Additional arguments to be passed to the generic methods.
  }
}
\value{
  Meta-community diversity objects are lists containing:
  \item{MetaCommunity}{The name of the \code{\link{MetaCommunity}} object containing inventory data.}
  \item{Type}{The type of diversity ("alpha", "beta" or "gamma").}
  \item{Order}{The order of diversity \eqn{q}.}
  \item{Correction}{The estimation bias correction used to calculate diversity.}
  \item{Tree}{The phylogenetic or functional tree used to calculate phylodiversity.}
  \item{Normalized}{Logical. Indicates whether phylodiversity is normalized or proportional to the height of the tree.}
  \item{Weights}{A vector containing the weights of communities.}
  \item{Communities}{A vector containing the diversity of communities.}
  \item{Total}{The total diversity.}

  \code{is.MCdiversity} returns \code{TRUE} if the object is of class \code{MCdiversity}.
  
  \code{summary.MCdiversity} returns a summary of the object's value.
}

\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
