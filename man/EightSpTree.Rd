\name{EightSpTree}
\docType{data}
\alias{EightSpTree}
\title{
  Functional tree with 8 species.
}
\description{
  This dataset is a leight-weight example.
}
\usage{
data(Paracou618)
}
\format{
  An object of class \code{\link[ade4]{phylog}} containing a functional tree.
}

\examples{
data(Paracou618)
# Preprocess the tree to be able to plot it 
# without loading ade4 package
plot(Preprocess.Tree(EightSpTree), hang=-0.01)
}

\keyword{datasets}
