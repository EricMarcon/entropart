# Tsallis beta entropy of a community

Calculates the Tsallis beta entropy of order \\q\\ of a community
belonging to a metacommunity.

## Usage

``` r
TsallisBeta(NorP, NorPexp = NULL, q = 1, ...)
bcTsallisBeta(Ns, Nexp = NULL, q, Correction = "Best", CheckArguments = TRUE)
# S3 method for class 'ProbaVector'
TsallisBeta(NorP, NorPexp = NULL, q = 1, ..., 
  CheckArguments = TRUE, Ps = NULL, Pexp = NULL)
# S3 method for class 'AbdVector'
TsallisBeta(NorP, NorPexp = NULL, q = 1, Correction = "Best", ..., 
  CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
# S3 method for class 'integer'
TsallisBeta(NorP, NorPexp = NULL, q = 1, Correction = "Best", ..., 
  CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
# S3 method for class 'numeric'
TsallisBeta(NorP, NorPexp = NULL, q = 1, Correction = "Best", ..., 
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

- q:

  A number: the order of entropy. Default is 1 for Shannon entropy.

- Correction:

  A string containing one of the possible corrections: currently, only
  `"ChaoShen"` or `"None"`. `"Best"` is the default value, it is
  equivalent to `"ChaoShen"`.

- ...:

  Additional arguments. Unused.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

The derivation of Tsallis beta entropy can be found in Marcon et al.
(2014).

Bias correction requires the number of individuals to estimate sample
[`Coverage`](https://ericmarcon.github.io/entropart/reference/Coverage.md).
Use `bcTsallisBeta` and choose the `Correction`.

Note that beta entropy value is related to alpha entropy (if \\q\\ is
not 1) and cannot be compared accross communities (Jost, 2007). Beta
entropy of a community is not meaningful in general, do rather calculate
the
[`BetaDiversity`](https://ericmarcon.github.io/entropart/reference/BetaDiversity.md)
of the metacommunity.

The functions are designed to be used as simply as possible.
`TsallisBeta` is a generic method. If its first argument is an abundance
vector, an integer vector or a numeric vector which does not sum to 1,
the bias corrected function `bcTsallisBeta` is called. Explicit calls to
`bcTsallisBeta` (with bias correction) or to `TsallisBeta.ProbaVector`
(without correction) are possible to avoid ambiguity. The `.integer` and
`.numeric` methods accept `Ps` or `Ns` arguments instead of `NorP` for
backward compatibility.

## Value

A named number equal to the calculated entropy. The name is that of the
bias correction used.

## References

Jost (2007), Partitioning diversity into independent alpha and beta
components. *Ecology* 88(10): 2427-2439.

Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014).
Generalization of the partitioning of Shannon diversity. *PLOS One*
9(3): e90289.

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- Paracou618.MC$Ps
# Probability distribution of the first plot
Ps1 <- Paracou618.MC$Psi[, 1]
# Divergence of order 2 between plot 1 and the whole forest
TsallisBeta(Ps1, Ps, 2)
#>       None 
#> 0.01004705 

# Ns is the vector of abundances of the metacommunity
Ns <- Paracou618.MC$Ns
# Abundances in the first plot
Ns1 <- Paracou618.MC$Nsi[, 1]
# Divergence of order 2 between plot 1 and the whole forest, with bias correction
bcTsallisBeta(Ns1, Ns, 2)
#>        Best 
#> 0.007566555 
```
