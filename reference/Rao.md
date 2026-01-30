# Rao Quadratic Entropy of a Community

Calculates Rao's quadratic entropy of a community described by a
probability vector and a phylogenetic / functional tree.

## Usage

``` r
Rao(NorP, Tree, ...)
bcRao(Ns, Tree, Correction="Lande", CheckArguments = TRUE)
# S3 method for class 'ProbaVector'
Rao(NorP, Tree, ..., CheckArguments = TRUE, Ps = NULL)
# S3 method for class 'AbdVector'
Rao(NorP, Tree, Correction = "Lande", ..., 
  CheckArguments = TRUE, Ns = NULL)
# S3 method for class 'integer'
Rao(NorP, Tree, Correction = "Lande", ..., 
  CheckArguments = TRUE, Ns = NULL)
# S3 method for class 'numeric'
Rao(NorP, Tree, Correction = "Lande", ..., 
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

- Tree:

  An object of class [`hclust`](https://rdrr.io/r/stats/hclust.html),
  "phylo" (see
  [`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)),
  [`phylog`](https://adeverse.github.io/ade4/reference/phylog.html) or
  [`PPtree`](https://ericmarcon.github.io/entropart/reference/PPtree.md).
  The tree must be ultrametric.

- Correction:

  A string containing one of the possible corrections accepted by
  [`bcTsallis`](https://ericmarcon.github.io/entropart/reference/Tsallis.md)
  or `"Lande"`, the default value (equivalent to `"Best"`).

- ...:

  Additional arguments. Unused.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

Bias correction requires the number of individuals. Use `bcRao` and
choose the `Correction`.

The unbiased estimator of Rao's entropy is identical to that of
Simpson's entropy because Rao's entropy is a linear sum of Simson
entropies, all of them calculated from the same number of individuals
(Marcon and Herault, 2014). It equals the plug-in etimator multiplied by
n/(n-1) where n is the total number of individuals.

The functions are designed to be used as simply as possible. `Tsallis`
is a generic method. If its first argument is an abundance vector, an
integer vector or a numeric vector which does not sum to 1, the bias
corrected function `bcTsallis` is called. Explicit calls to `bcTsallis`
(with bias correction) or to `Tsallis.ProbaVector` (without correction)
are possible to avoid ambiguity. The `.integer` and `.numeric` methods
accept `Ps` or `Ns` arguments instead of `NorP` for backward
compatibility.

## Value

A named number equal to the calculated entropy. The name is that of the
bias correction used.

## References

Marcon, E., Herault, B. (2015). Decomposing Phylodiversity. *Methods in
Ecology and Evolution* 6(3): 333-339.

Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified
approach. *Theoretical Population Biology* 21: 24-43.

## See also

[`bcPhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Species probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ns)
# Calculate Rao's quadratic entropy of the plot
Rao(Ps, Paracou618.Taxonomy)
#>     None 
#> 2.878133 
```
