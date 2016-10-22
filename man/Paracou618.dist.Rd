\name{Paracou618.dist}
\docType{data}
\alias{Paracou618.dist}
\title{
  Functional distances between pairs of species of Paracou field station plots 6 and 18, two 1-ha plots inventoried by the Bridge project.
}
\description{
  This dataset is from Paracou field station, French Guiana, managed by \href{http://www.cirad.fr}{Cirad}. Traits are detailed in Marcon and Herault (2014), the distance matrix was built following Paine et al. (2011).
}
\usage{
data(Paracou618)
}
\format{
  An object of class \code{\link{dist}}.
}

\source{
  Permanent data census of Paracou.
}
\references{
  Gourlet-Fleury, S., Guehl, J. M. and Laroussinie, O., Eds. (2004). \emph{Ecology & management of a neotropical rainforest. Lessons drawn from Paracou, a long-term experimental research site in French Guiana}. Paris, Elsevier.
  
  Marcon, E., Herault, B. (2015). Decomposing  Phylodiversity. \emph{Methods in Ecology and Evolution} 6(3): 333-339.
  
  Paine, C. E. T., Baraloto, C., Chave, J., and Herault, B. (2011). Functional traits of individual trees reveal ecological constraints on community assembly in tropical rain forests. \emph{Oikos}, 120(5), 720-727.
}
\examples{
data(Paracou618)
plot(density(Paracou618.dist, from=0), main="Distances between species")
}

\keyword{datasets}
