# Similarity-based diversity of a community

Calculates the diversity of order \\q\\ of a probability vector
according to a similarity matrix.

## Usage

``` r
Dqz(NorP, q = 1, Z = diag(length(NorP)), ...)
bcDqz(Ns, q = 1, Z = diag(length(Ns)), Correction = "Best", CheckArguments = TRUE)
# S3 method for class 'ProbaVector'
Dqz(NorP, q = 1, Z = diag(length(NorP)), ..., 
  CheckArguments = TRUE, Ps = NULL)
# S3 method for class 'AbdVector'
Dqz(NorP, q = 1, Z = diag(length(NorP)), Correction = "Best", ..., 
  CheckArguments = TRUE, Ns = NULL)
# S3 method for class 'integer'
Dqz(NorP, q = 1, Z = diag(length(NorP)), Correction = "Best", ..., 
  CheckArguments = TRUE, Ns = NULL)
# S3 method for class 'numeric'
Dqz(NorP, q = 1, Z = diag(length(NorP)), Correction = "Best", ..., 
  CheckArguments = TRUE, Ps = NULL, Ns = NULL)
```

## Arguments

- Ps:

  A probability vector, summing to 1.

- Ns:

  A numeric vector containing species abundances.

- NorP:

  A numeric vector, an integer vector, an abundance vector
  ([`AbdVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md))
  or a probability vector
  ([`ProbaVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md)).
  Contains either abundances or probabilities.

- q:

  A number: the order of diversity. Default is 1.

- Z:

  A relatedness matrix, *i.e.* a square matrix whose terms are all
  positive, strictly positive on the diagonal. Generally, the matrix is
  a similarity matrix, *i.e.* the diagonal terms equal 1 and other terms
  are between 0 and 1. Default is the identity matrix to calculate
  neutral diversity.

- Correction:

  A string containing one of the possible corrections: `"None"` (no
  correction), `"HorvitzThomson"`, `"MarconZhang"` or `"Best"`, the
  default value. The `"MarconZhang"` correction assumes a similarity
  matrix.

- ...:

  Additional arguments. Unused.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

Diversity is calculated following Leinster and Cobbold (2012): it is the
reciprocal of the (generalized) average (of order `q`) of the community
species ordinariness.

A similarity matrix is used (as for `Dqz`), not a distance matrix as in
Ricotta and Szeidl (2006). See the example.

Bias correction requires the number of individuals. Use `bcHqz` and
choose the `Correction`. Correction techniques are from Marcon *et al.*
(2014).

Currently, the `"Best"` correction is the max value of
`"HorvitzThomson"` and `"MarconZhang"`.

The functions are designed to be used as simply as possible. `Dqz` is a
generic method. If its first argument is an abundance vector, an integer
vector or a numeric vector which does not sum to 1, the bias corrected
function `bcDqz` is called. Explicit calls to `bcDqz` (with bias
correction) or to `Dqz.ProbaVector` (without correction) are possible to
avoid ambiguity. The `.integer` and `.numeric` methods accept `Ps` or
`Ns` arguments instead of `NorP` for backward compatibility.

## Value

A named number equal to the calculated diversity. The name is that of
the bias correction used.

## References

Leinster, T. and Cobbold, C. (2012). Measuring diversity: the importance
of species similarity. *Ecology* 93(3): 477-489.

Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of
similarity-based diversity and its bias correction. *HAL*
hal-00989454(version 3).

## See also

[`Hqz`](https://ericmarcon.github.io/entropart/reference/Hqz.md),
[`PhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Prepare the similarity matrix
DistanceMatrix <- as.matrix(Paracou618.dist)
# Similarity can be 1 minus normalized distances between species
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
# Calculate diversity of order 2
Dqz(Paracou618.MC$Ns, 2, Z)
#>    Best 
#> 1.48295 
```
