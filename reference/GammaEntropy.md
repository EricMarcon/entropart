# Reduced-bias gamma entropy of a metacommunity

Calculates the reduced-bias Tsallis entropy of order \\q\\ of a
metacommunity.

## Usage

``` r
GammaEntropy(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, 
  Z = NULL, PhyloDetails = FALSE, CheckArguments = TRUE)
```

## Arguments

- MC:

  A
  [`MetaCommunity`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md)
  object.

- q:

  A number: the order of entropy. Default is 1.

- Correction:

  A string containing one of the possible corrections accepted by the
  bias-corrected entropy function (see details) or `"None"` or `"Best"`,
  the default value.

- Tree:

  An object of class [`hclust`](https://rdrr.io/r/stats/hclust.html),
  "phylo" (see
  [`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)),
  [`phylog`](https://adeverse.github.io/ade4/reference/phylog.html) or
  [`PPtree`](https://ericmarcon.github.io/entropart/reference/PPtree.md).
  The tree must be ultrametric.

- Normalize:

  If `TRUE` (default), the entropy returned by the function is
  normalized by the height of the tree (it is the weighted average value
  of the entropy in each slice).  
  If `FALSE`, it is the unnormalized weighted sum of the results.

- Z:

  A relatedness matrix, *i.e.* a square matrix whose terms are all
  positive, strictly positive on the diagonal. Generally, the matrix is
  a similarity matrix, *i.e.* the diagonal terms equal 1 and other terms
  are between 0 and 1.

- PhyloDetails:

  If `FALSE` (default), the function always returns a number. If `TRUE`
  and `Tree` is not Tree is not `NULL` then a
  [`PhyloValue`](https://ericmarcon.github.io/entropart/reference/PhyloValue.md)
  object is returned with all details. That is used internally by
  [`DivPart`](https://ericmarcon.github.io/entropart/reference/DivPart.md)
  to obtain the corrections used to estimate gamma entropy along the
  tree and apply them to the estimation of alpha diversity.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

If `Tree` is not `NULL`, then phylogenetic entropy is calculated by
[`bcPhyloEntropy`](https://ericmarcon.github.io/entropart/reference/PhyloEntropy.md).

Else, if `Z` is not `NULL`, then similarity-based entropy is calculated
by [`bcHqz`](https://ericmarcon.github.io/entropart/reference/Hqz.md).

Else, neutral entropy is calculated by
[`bcTsallis`](https://ericmarcon.github.io/entropart/reference/Tsallis.md).

## Value

A number equal to the calculated entropy.

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

[`bcTsallis`](https://ericmarcon.github.io/entropart/reference/Tsallis.md),
[`bcPhyloEntropy`](https://ericmarcon.github.io/entropart/reference/PhyloEntropy.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Calculate Simpson gamma entropy
GammaEntropy(Paracou618.MC, 2)
#>   UnveilJ 
#> 0.9863189 
# Compare without correction
GammaEntropy(Paracou618.MC, 2, Correction = "None")
#>      None 
#> 0.9854485 
# Estimate phylogenetic Simpson gamma entropy
GammaEntropy(Paracou618.MC, 2, Tree = Paracou618.Taxonomy)   
#> [1] 0.9601824
```
