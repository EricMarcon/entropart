\name{ArgumentOriginalName}
\alias{ArgumentOriginalName}
\title{
  Finds the original value (name of expression) of the argument of a function even in the case of embedded calls.
}
\description{
  This function is used internally to find the name of arguments passed to entropart functions such as \code{\link{PhyloDiversity}} that store them in their results.
}
\usage{
ArgumentOriginalName(x)
}
\arguments{
  \item{x}{
  Any argument.
  }
}
\value{
  The name of the argument.
}
\details{
  The function searches the name of the argument in the parent frame of the function and stops at the top.
}
\author{
  BrodieG, <https://stackoverflow.com/users/2725969/brodieg>
  
  Eric Marcon <Eric.Marcon@agroparistech.fr>
}


\keyword{internal}
