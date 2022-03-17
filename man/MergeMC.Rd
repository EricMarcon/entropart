\name{MC Utilities}
\alias{MergeMC}
\alias{MergeC}
\alias{ShuffleMC}
\title{
  Manipulation of meta-communities
}
\description{
  Tools to manipulate meta-communities. From a list of meta-communities, \code{MergeMC} creates a metacommunity whose communities are each original metacommunity. \code{MergeC} creates a metacommunity whose communities are each original community. \code{ShuffleMC} randomly assigns original communities to a metacommunity, keeping original weights, and returns a list of meta-communities.
}
\usage{
MergeMC(MClist, Weights = rep(1, length(MClist)), CheckArguments = TRUE)
MergeC(MClist, Weights = rep(1, length(MClist)), CheckArguments = TRUE)
ShuffleMC(MClist, Weights = rep(1, length(MClist)), CheckArguments = TRUE)
}
\arguments{
  \item{MClist}{
  A list of \code{MetaCommunity} objects.
  }
  \item{Weights}{
  A vector of numbers containing the weight of each metacommunity of the list. It does not have to be normalized to sum to 1.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  \code{MergeMC} is used for hierarchical partitioning of diversity. The gamma diversity of communities of the list becomes alpha diversity of the merged meta-community.
  
  \code{MergeC} creates a new meta-community by mixing original ones. Original communities are kept, their new weight is the product of their original weight and the weight of their original meta-community.
  
  \code{ShuffleMC} is used for simulations of the null hypothesis that all metacommunities of the list are identical.
}
\value{
  \code{MergeMC} and \code{MergeC} return a \code{\link{MetaCommunity}}.
  
  \code{ShuffleMC} returns a \code{\link{list}} of \code{\link{MetaCommunity}} objects.
}
\seealso{
  \code{\link{MetaCommunity}}
}
\examples{
# First meta-community
(df <- data.frame(C1 = c(10, 10, 10, 10), C2 = c(0, 20, 35, 5),
  C3 = c(25, 15, 0, 2), row.names = c("sp1", "sp2", "sp3", "sp4")))
w <- c(1, 2, 1)
MC1 <- MetaCommunity(Abundances = df, Weights = w)
# Second meta-community
(df <- data.frame(C1 = c(10, 4), C2 = c(3, 4), row.names = c("sp1", "sp5")))
w <- c(3, 2)
MC2 <- MetaCommunity(Abundances = df, Weights = w)

# Merge communities
plot(MergeC(list(MC1, MC2)), main="Merged communities")
# Merge metacommunities
plot(MergeMC(list(MC1, MC2)), main="Merged meta-communities")
smc <- ShuffleMC(list(MC1, MC2))
plot(MergeMC(smc), main="Shuffled, then Merged meta-communities")
}
