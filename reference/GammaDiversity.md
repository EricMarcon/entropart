# Reduced-bias gamma diversity of a metacommunity

Calculates the reduced-bias diversity of order \\q\\ of a metacommunity.

## Usage

``` r
GammaDiversity(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, 
  Z = NULL, CheckArguments = TRUE)
```

## Arguments

- MC:

  A
  [`MetaCommunity`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md)
  object.

- q:

  A number: the order of diversity. Default is 1.

- Correction:

  A string containing one of the possible corrections accepted by
  [`AlphaEntropy`](https://ericmarcon.github.io/entropart/reference/AlphaEntropy.md)
  or `"None"` or `"Best"`, the default value.

- Tree:

  An object of class [`hclust`](https://rdrr.io/r/stats/hclust.html),
  "phylo" (see
  [`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)),
  [`phylog`](https://adeverse.github.io/ade4/reference/phylog.html) or
  [`PPtree`](https://ericmarcon.github.io/entropart/reference/PPtree.md).
  The tree must be ultrametric.

- Normalize:

  If `TRUE` (default), diversity is not affected by the height of the
  tree.  
  If `FALSE`, diversity is proportional to the height of the tree.

- Z:

  A relatedness matrix, *i.e.* a square matrix whose terms are all
  positive, strictly positive on the diagonal. Generally, the matrix is
  a similarity matrix, *i.e.* the diagonal terms equal 1 and other terms
  are between 0 and 1.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

Entropy is calculated by
[`GammaEntropy`](https://ericmarcon.github.io/entropart/reference/GammaEntropy.md)
and transformed into diversity.

## Value

The metacommunity's gamma entropy.

## References

Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014).
Generalization of the partitioning of Shannon diversity. *PLOS One*
9(3): e90289.

Marcon, E., Herault, B. (2015). Decomposing Phylodiversity. *Methods in
Ecology and Evolution* 6(3): 333-339.

Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of
similarity-based diversity and its bias correction. *HAL*
hal-00989454(version 3).

## See also

[`GammaEntropy`](https://ericmarcon.github.io/entropart/reference/GammaEntropy.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Calculate Simpson gamma diversity
GammaDiversity(Paracou618.MC, 2)
#>  UnveilJ 
#> 73.09335 
# Compare without correction
GammaDiversity(Paracou618.MC, 2, Correction = "None")
#>    None 
#> 68.7215 
# Estimate phylogenetic Simpson gamma diversity
GammaDiversity(Paracou618.MC, 2, Tree = Paracou618.Taxonomy)   
#> [1] 25.11451
```
