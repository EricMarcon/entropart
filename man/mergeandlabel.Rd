\name{mergeandlabel}
\alias{mergeandlabel}
\title{
  Merge community data
}
\description{
  Merge two dataframes containing species abundances.
}
\usage{
mergeandlabel(x, y)
}
\arguments{
  \item{x}{
  A dataframe containing species abundances
  }
  \item{y}{
  A dataframe containing species abundances
  }
}
\details{
  Row names must contain species. The new dataframes keeps all species with abundance zero for species that were not found in orfginal communities. 
}
\value{
  A dataframe containing all lines and columns of the merged dataframes.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{MergeC}}, \code{\link{MergeMC}}, \code{\link{ShuffleMC}}
}
\keyword{internal}
