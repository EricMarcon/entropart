# Similarity-based beta entropy of a community

Calculates the similarity-based beta entropy of order \\q\\ of a
community belonging to a metacommunity.

## Usage

``` r
HqzBeta(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), ...)
bcHqzBeta(Ns, Nexp = NULL, q = 1, Z = diag(length(Ns)), Correction = "Best",
          CheckArguments = TRUE)
# S3 method for class 'ProbaVector'
HqzBeta(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)),
  ..., CheckArguments = TRUE, Ps = NULL, Pexp = NULL)
# S3 method for class 'AbdVector'
HqzBeta(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), Correction = "Best", 
  ..., CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
# S3 method for class 'integer'
HqzBeta(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), Correction = "Best", 
  ..., CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
# S3 method for class 'numeric'
HqzBeta(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), Correction = "Best", 
  ..., CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL)
```

## Arguments

- Ps:

  The probability vector of species of the community.

- Pexp:

  The probability vector of species of the metacommunity.

- Ns:

  A numeric vector containing species abundances of the community.

- Nexp:

  A numeric vector containing species abundances of the metacommunity.

- NorP:

  A numeric vector, an integer vector, an abundance vector
  ([`AbdVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md))
  or a probability vector
  ([`ProbaVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md)).
  Contains either abundances or probabilities of the community.

- NorPexp:

  A numeric vector, an integer vector, an abundance vector
  ([`AbdVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md))
  or a probability vector
  ([`ProbaVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md)).
  Contains either abundances or probabilities of the metacommunity.

- q:

  A number, the order of diversity. Default is 1.

- Z:

  A relatedness matrix, *i.e.* a square matrix whose terms are all
  positive, strictly positive on the diagonal. Generally, the matrix is
  a similarity matrix, *i.e.* the diagonal terms equal 1 and other terms
  are between 0 and 1. Default is the identity matrix to calculate
  neutral entropy.

- Correction:

  A string containing one of the possible corrections: currently, no
  correction is available so `"Best"`, the default value, is equivalent
  to `"None"`.

- ...:

  Additional arguments. Unused.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

The derivation of similarity-based beta entropy can be found in Marcon
*et al.* (2014).

Bias correction requires the number of individuals.

Note that beta entropy value is related to alpha entropy (if \\q\\ is
not 1) and cannot be compared accross communities (Jost, 2007). Beta
entropy of a community is not meaningful in general, do rather calculate
the
[`BetaDiversity`](https://ericmarcon.github.io/entropart/reference/BetaDiversity.md)
of the metacommunity.

The functions are designed to be used as simply as possible. `HqzBeta`
is a generic method. If its first argument is an abundance vector, an
integer vector or a numeric vector which does not sum to 1, the bias
corrected function `bcHqzBeta` is called. Explicit calls to `bcHqzBeta`
(with bias correction) or to `HqzBeta.ProbaVector` (without correction)
are possible to avoid ambiguity. The `.integer` and `.numeric` methods
accept `Ps` or `Ns` arguments instead of `NorP` for backward
compatibility.

## Value

A named number equal to the calculated entropy. The name is that of the
bias correction used.

## References

Jost (2007), Partitioning diversity into independent alpha and beta
components. *Ecology* 88(10): 2427-2439.

Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of
similarity-based diversity and its bias correction. *HAL*
hal-00989454(version 3).

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ps)
# Probability distribution of the first plot
Ps1 <- as.ProbaVector(Paracou618.MC$Psi[, 1])
# Prepare the similarity matrix
DistanceMatrix <- as.matrix(Paracou618.dist)
# Similarity can be 1 minus normalized distances between species
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
# Divergence of order 2 between plot 1 and the whole forest
HqzBeta(Ps1, Ps, q=2, Z)
#>        None 
#> 0.004107207 
```
