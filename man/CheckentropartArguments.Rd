\name{CheckentropartArguments}
\alias{CheckentropartArguments}
\title{
  Checks the arguments of a function of the package entropart
}
\description{
  This function is used internally to verify that arguments passed to entropart functions such as \code{\link{PhyloDiversity}} are correct.
}
\usage{
CheckentropartArguments()
}
\value{
  Returns \code{TRUE} or stops if a problem is detected.
}
\details{
  The function compares the arguments passed to its parent function to the type they should be and performs some extra tests (\emph{e.g.} probabilities must be positive and sum to 1). It stops if an argument is not correct.
}
\keyword{internal}
