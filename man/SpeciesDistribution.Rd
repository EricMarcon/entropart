\name{SpeciesDistribution}
\alias{SpeciesDistribution}
\alias{as.SpeciesDistribution}
\alias{as.SpeciesDistribution.data.frame}
\alias{as.SpeciesDistribution.integer}
\alias{as.SpeciesDistribution.numeric}
\alias{is.SpeciesDistribution}
\alias{plot.SpeciesDistribution}
\alias{AbdVector}
\alias{as.AbdVector}
\alias{as.AbdVector.data.frame}
\alias{as.AbdVector.integer}
\alias{as.AbdVector.numeric}
\alias{is.AbdVector}
\alias{ProbaVector}
\alias{as.ProbaVector}
\alias{as.ProbaVector.data.frame}
\alias{as.ProbaVector.integer}
\alias{as.ProbaVector.numeric}
\alias{is.ProbaVector}
\title{
  Species Distributions.
}
\description{
  A Species Distribution is a (preferably named) vector containing species abundances or probabilities.
}
\usage{
as.SpeciesDistribution(x)
\method{as.SpeciesDistribution}{data.frame}(x)
\method{as.SpeciesDistribution}{integer}(x)
\method{as.SpeciesDistribution}{numeric}(x)
\method{plot}{SpeciesDistribution}(x, \dots, Distribution = NULL, 
         type = "b", log = "y", main = NULL, xlab = "Rank", ylab = NULL)
is.SpeciesDistribution(x)
as.ProbaVector(x, Correction = "None", Unveiling = "None", 
          RCorrection = "Chao1", JackOver = FALSE, CEstimator = "ZhangHuang", 
          CheckArguments = TRUE)
\method{as.ProbaVector}{data.frame}(x, Correction = "None", Unveiling = "None", 
          RCorrection = "Chao1", JackOver = FALSE, CEstimator = "ZhangHuang", 
          CheckArguments = TRUE)
\method{as.ProbaVector}{integer}(x, Correction = "None", Unveiling = "None", 
          RCorrection = "Chao1", JackOver = FALSE, CEstimator = "ZhangHuang", 
          CheckArguments = TRUE)
\method{as.ProbaVector}{numeric}(x, Correction = "None", Unveiling = "None", 
          RCorrection = "Chao1", JackOver = FALSE, CEstimator = "ZhangHuang", 
          CheckArguments = TRUE)
is.ProbaVector(x)
as.AbdVector(x, Round = TRUE)
\method{as.AbdVector}{data.frame}(x, Round = TRUE)
\method{as.AbdVector}{integer}(x, Round = TRUE)
\method{as.AbdVector}{numeric}(x, Round = TRUE)
is.AbdVector(x)
}
\arguments{
  \item{x}{
    An object.
  }
  \item{Distribution}{
    The distribution to fit on the plot. May be \code{"lnorm"} (log-normal), \code{"lseries"} (log-series), \code{"geom"} (geometric) or \code{"bstick"} (broken stick). If \code{NULL}, no distribution is fitted. See \code{\link{rCommunity}} for the description of these distributions.
  }
  \item{Round}{
    If \code{TRUE} (by default), values of \code{x} are set to integer to create an \code{AbdVector}. This is useful if original abundances are not integers (this is often the case for \code{\link{MetaCommunity}} abundances which are the product of probabilities by the number of individuals) and integer values are required (for example to calculate the bootstrap confidence interval of a community profile).
  }
  \item{Correction}{
    A string containing one of the possible corrections to estimate a probability distribution: \code{"None"} (no correction, the default value), \code{"Chao2013"}, \code{"Chao2015"} or \code{"ChaoShen"}.
  }
  \item{Unveiling}{
    A string containing one of the possible unveiling methods to estimate the probabilities of the unobserved species: \code{"None"} (default, no species is added), \code{"unif"} (uniform: all unobserved species have the same probability) or \code{"geom"} (geometric: the unobserved species distribution is geometric).
  }
  \item{RCorrection}{
    A string containing a correction recognized by \code{\link{Richness}} to evaluate the total number of species. \code{"Chao1"} is  the default value.
  }
  \item{JackOver}{
    If \code{TRUE}, retain the jackknife order immediately superior to the optimal one, usually resulting in the overestimation of the number of species. Default is \code{FALSE}. Ignored if \code{RCorrection} is not \code{"Jackknife"}.
  }
  \item{CEstimator}{
    A string containing an estimator recognized by \code{\link{Coverage}} to evaluate the sample coverage. \code{"ZhangHuang"} is the default value.
  }
 \item{type}{
    The plot type, see \code{\link{plot}}.
  }
  \item{log}{
    The axis to plot in log scale, \emph{e.g.} \code{"xy"} for both axes. Default is \code{"y"}.
  }
  \item{main}{
    The main title of the plot. if \code{NULL} (by default), there is no title.
  }
  \item{xlab}{
    The X axis label, "Rank" by default.
  }
  \item{ylab}{
    The Y axis label. if \code{NULL} (by default), "Probability" or "Abundance" is chosen according to the object class.
  }
  \item{\dots}{
    Additional arguments to be passed to \code{\link{plot}}.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  \code{SpeciesDistribution} objects include \code{AbdVector} and \code{ProbaVector} objects.
  
  \code{as.AbdVector} just sets the class of the numeric or integer \code{x} so that appropriate versions of community functions (generic methods such as \code{\link{Diversity}}) are applied. Abundance values are rounded (by default) to the nearest integer.
  
  \code{as.ProbaVector} normalizes the vector so that it sums to 1. If \code{Correction} is not \code{"None"}, the observed abundance distribution is used to estimate the actual species distribution. The list of species will be changed: zero-aundance species will be cleared, and some unobserved species will be added. First, observed species probabilities are estimated folllowing Chao and Shen (2003), \emph{i.e.} input probabilities are multiplied by the sample coverage, or according to more sophisticated models: Chao \emph{et al.} (2013, single-parameter model), or Chao \emph{et al.} (2015, two-parameter model). The total probability of observed species equals the sample coverage. Then, the distribution of unobserved species can be unveiled: their number is estimated according to \code{RCorrection} (if the Jackknife estimator is chosen, the \code{JackOver} argument allows using the order immediately over the optimal one). The coverage deficit (1 minus the sample coverage) is shared by the unobserved species equally (\code{Unveiling = "unif"}, Chao \emph{et al.}, 2013) or according to a geometric distribution (\code{Unveiling = "geom"}, Chao \emph{et al.}, 2015).
  
  These functions can be applied to data frames to calculate the joint diversity (Gregorius, 2010).
  
  \code{SpeciesDistribution} objects can be plotted. The \code{plot} method returns the estimated parameters of the fitted distribution. The broken stick has no parameter, so the maximum abundance is returned.
}
\note{
  Fisher's alpha (Fisher \emph{et al.}, 1943) is estimated to fit the log-series distribution. The estimation is done by the \code{\link{fisher.alpha}} function of package \code{vegan}. It may differ substantially from the estimation returned by \code{optimal.theta} from package \code{untb}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>, Bruno Herault <Bruno.Herault@ecofog.gf>
}
\references{
  Chao, A. and Shen, T. J. (2003). Nonparametric estimation of Shannon's index of diversity when there are unseen species in sample. \emph{Environmental and Ecological Statistics} 10(4): 429-443.

  Chao, A., Wang, Y. T. and Jost, L. (2013). Entropy and the species accumulation curve: a novel entropy estimator via discovery rates of new species. \emph{Methods in Ecology and Evolution} 4(11):1091-1100.

  Chao, A., Hsieh, T. C., Chazdon, R. L., Colwell, R. K., Gotelli, N. J. (2015) Unveiling the Species-Rank Abundance Distribution by Generalizing Good-Turing Sample Coverage Theory. \emph{Ecology} 96(5): 1189-1201.

  Fisher R.A., Corbet A.S., Williams C.B. (1943) The Relation Between the Number of Species and the Number of Individuals in a Random Sample of an Animal Population. \emph{Journal of Animal Ecology} 12: 42-58.

  Gregorius H.-R. (2010) Linking Diversity and Differentiation. \emph{Diversity} 2(3): 370-394.
}
\seealso{
  \code{\link{rgeom}}, \code{\link{rlnorm}}, \code{\link{rCommunity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ns is the total number of trees per species
  Ns <- as.AbdVector(Paracou618.MC$Ns)
  # Whittaker plot, poorly fitted by a log-normal distribution
  plot(Ns, Distribution = "lnorm")
}
