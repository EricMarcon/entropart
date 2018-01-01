\name{Paracou618.Taxonomy}
\docType{data}
\alias{Paracou618.Taxonomy}
\title{
  Taxonomy (Family - Genus - Species) of Paracou field station plots 6 and 18, two 1-ha plots inventoried by the Bridge project.
}
\description{
  This dataset is from Paracou field station, French Guiana, managed by \href{http://www.cirad.fr}{Cirad}.
}
\usage{
data(Paracou618)
}
\format{
  An object of class \code{\link{phylo}} containing a taxonomy.
}

\source{
  Permanent data census of Paracou.
}
\references{
  Gourlet-Fleury, S., Guehl, J. M. and Laroussinie, O., Eds. (2004). \emph{Ecology & management of a neotropical rainforest. Lessons drawn from Paracou, a long-term experimental research site in French Guiana}. Paris, Elsevier.
}
\examples{
data(Paracou618)
plot(Paracou618.Taxonomy, type="fan", show.tip.label=FALSE)
}

\keyword{datasets}
