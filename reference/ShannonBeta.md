# Shannon beta entropy of a community

Calculates the Shannon beta entropy of a community belonging to a
metacommunity.

## Usage

``` r
ShannonBeta(NorP, NorPexp = NULL, ...)
bcShannonBeta(Ns, Nexp, Correction = "Best", CheckArguments = TRUE)
# S3 method for class 'ProbaVector'
ShannonBeta(NorP, NorPexp = NULL, ..., 
  CheckArguments = TRUE, Ps = NULL, Pexp = NULL)
# S3 method for class 'AbdVector'
ShannonBeta(NorP, NorPexp = NULL, Correction = "Best", ...,
  CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
# S3 method for class 'integer'
ShannonBeta(NorP, NorPexp = NULL, Correction = "Best", ...,
  CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
# S3 method for class 'numeric'
ShannonBeta(NorP, NorPexp = NULL, Correction = "Best", ...,
  CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL)
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

- Correction:

  A string containing one of the possible corrections: currently,
  `"ChaoShen"` (Marcon *et al.*, 2012) equivalent to `"Best"`, and
  `"ZhangGrabchak"` (Zhang and Grabchak, 2014).

- ...:

  Additional arguments. Unused.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

The derivation of Shannon beta entropy can be found in Marcon *et al.*
(2012).

Bias correction requires the number of individuals to estimate sample
[`Coverage`](https://ericmarcon.github.io/entropart/reference/Coverage.md).
Use `bcShannonBeta` and choose the `Correction`.

The functions are designed to be used as simply as possible.
`ShannonBeta` is a generic method. If its first argument is an abundance
vector, an integer vector or a numeric vector which does not sum to 1,
the bias corrected function `bcShannonBeta` is called. Explicit calls to
`bcShannonBeta` (with bias correction) or to `ShannonBeta.ProbaVector`
(without correction) are possible to avoid ambiguity. The `.integer` and
`.numeric` methods accept `Ps` or `Ns` arguments instead of `NorP` for
backward compatibility.

## Value

A number equal to the calculated entropy.

## References

Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The
Decomposition of Shannon's Entropy and a Confidence Interval for Beta
Diversity. *Oikos* 121(4): 516-522.

Zhang, Z. and Grabchak M. (2014). Nonparametric Estimation of
Kullback-Leibler Divergence. *Neural computation* 26(11): 2570-2593.

## See also

`bcShannonBeta`

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ps)
# Probability distribution of the first plot
Ps1 <- as.ProbaVector(Paracou618.MC$Psi[, 1])
# Shannon beta entropy of the plot
ShannonBeta(Ps1, Ps)
#>      None 
#> 0.3499358 

# Ns is the vector of abundances of the metacommunity
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Abundances in the first plot
Ns1 <- as.AbdVector(Paracou618.MC$Nsi[, 1])
# Reduced-bias estimator of Shannon beta entropy of the plot
bcShannonBeta(Ns1, Ns)
#>      Best 
#> 0.2949188 
```
