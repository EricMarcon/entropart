# Simpson beta entropy of a community

Calculates the Simpson beta entropy of a community belonging to a
metacommunity.

## Usage

``` r
SimpsonBeta(NorP, NorPexp = NULL, ...)
bcSimpsonBeta(Ns, Nexp, Correction = "Best", CheckArguments = TRUE)
# S3 method for class 'ProbaVector'
SimpsonBeta(NorP, NorPexp = NULL, ..., 
  CheckArguments = TRUE, Ps = NULL, Pexp = NULL)
# S3 method for class 'AbdVector'
SimpsonBeta(NorP, NorPexp = NULL, Correction = "Best", ..., 
  CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
# S3 method for class 'integer'
SimpsonBeta(NorP, NorPexp = NULL, Correction = "Best", ..., 
  CheckArguments = TRUE, Ns = NULL, Nexp = NULL)
# S3 method for class 'numeric'
SimpsonBeta(NorP, NorPexp = NULL, Correction = "Best", ..., 
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

  A string containing one of the possible corrections: currently, only
  `"ChaoShen"`, identical to `"Best"`.

- ...:

  Additional arguments. Unused.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

The derivation of Tsallis beta entropy (Simpson is Tsallis of order 2)
can be found in Marcon et al. (2014).

Bias correction requires the number of individuals to estimate sample
[`Coverage`](https://ericmarcon.github.io/entropart/reference/Coverage.md).
Use `bcSimpsonBeta` and choose the `Correction`.

Note that Simpson beta entropy value is related to Simpson alpha entropy
value and cannot be compared accross communities (Jost, 2007). Beta
entropy of a community is not meaningful in general, do rather calculate
the
[`BetaDiversity`](https://ericmarcon.github.io/entropart/reference/BetaDiversity.md)
of order 2 of the metacommunity.

The functions are designed to be used as simply as possible.
`SimpsonBeta` is a generic method. If its first argument is an abundance
vector, an integer vector or a numeric vector which does not sum to 1,
the bias corrected function `bcSimpsonBeta` is called. Explicit calls to
`bcSimpsonBeta` (with bias correction) or to `SimpsonBeta.ProbaVector`
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

## See also

[`Simpson`](https://ericmarcon.github.io/entropart/reference/Simpson.md),
`bcSimpsonBeta`,
[`BetaDiversity`](https://ericmarcon.github.io/entropart/reference/BetaDiversity.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ps)
# Probability distribution of the first plot
Ps1 <- as.ProbaVector(Paracou618.MC$Psi[, 1])

# Simpson beta entropy of the plot
SimpsonBeta(Ps1, Ps)
#>       None 
#> 0.01004705 
# Transform into diversity
expq(SimpsonBeta(Ps1, Ps)/(1-Simpson(Ps1)), 2)
#>     None 
#> 1.583005 

# Ns is the vector of abundances of the metacommunity
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Abundances in the first plot
Ns1 <- as.AbdVector(Paracou618.MC$Nsi[, 1])
# Reduced-bias Shannon beta entropy of the plot
bcSimpsonBeta(Ns1, Ns)
#>        Best 
#> 0.007566555 
```
