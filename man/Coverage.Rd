\name{Coverage}
\alias{Coverage}
\title{
  Sample coverage of a community
}
\description{
  Calculates an estimator of the sample coverage of a community described by its abundance vector.
}
\usage{
Coverage(Ns, Estimator = "ZhangHuang", CheckArguments = TRUE)
}
\arguments{
  \item{Ns}{
  A numeric vector containing species abundances.
  }
  \item{Estimator}{
  A string containing one of the possible estimators: \code{"ZhangHuang"} (default), \code{"Chao"}, \code{"Turing"}.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The sample coverage \eqn{C} of a community is the total probability of occurence of the species observed in the sample. \eqn{1-C} is the probability for an individual of the whole community to belong to a species that has not been sampled.
  
  The historical estimator is due to Turing (Good, 1953). It only relies on singletons (species observed only once). Chao's (Chao and Shen, 2010) estimator uses doubletons too and Zhang-Huang's (Chao et al., 1988; Zhang and Huang, 2007) uses the whole distribution.
}
\value{
  A named number equal to the calculated sample coverage. The name is that of the estimator used.
}
\references{
  Chao, A., Lee, S.-M. and Chen, T.-C. (1988). A generalized Good's nonparametric coverage estimator. \emph{Chinese Journal of Mathematics} 16: 189-199.
  
  Chao, A. and Shen, T.-J. (2010). \emph{Program SPADE: Species Prediction And Diversity Estimation. Program and user's guide}. CARE, Hsin-Chu, Taiwan.
  
  Good, I. J. (1953). On the Population Frequency of Species and the Estimation of Population Parameters. \emph{Biometrika} 40(3/4): 237-264.
  
  Zhang, Z. and Huang, H. (2007). Turing's formula revisited. \emph{Journal of Quantitative Linguistics} 14(2-3): 222-241.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the vector of abundances of the metacommunity
Ns <- Paracou618.MC$Ns
# Calculate the sample coverage of the metacommunity
Coverage(Ns)    # Stored in Paracou618.SampleCoverage
}
