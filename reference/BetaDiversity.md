# Reduced-bias beta diversity of a metacommunity

Calculates the reduced-bias beta diversity of order \\q\\ between
communities.

## Usage

``` r
BetaDiversity(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, 
  Z = NULL, CheckArguments = TRUE)
```

## Arguments

- MC:

  A
  [`MetaCommunity`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md)
  object.

- q:

  A number: the order of diversity. Default is 1 for Shannon diversity.

- Correction:

  A string containing one of the possible corrections accepted by
  [`bcTsallisBeta`](https://ericmarcon.github.io/entropart/reference/TsallisBeta.md)
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
[`BetaEntropy`](https://ericmarcon.github.io/entropart/reference/BetaEntropy.md)
and transformed into diversity.

Diversity values of communities are not defined: community entropies are
averaged to obtain the metacommunity entropy wich is transformed into
diversity (Marcon et al., 2014).

## Value

An
[`MCdiversity`](https://ericmarcon.github.io/entropart/reference/MCdiversity.md)
object containing diversity value of the metacommunity.

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

[`BetaEntropy`](https://ericmarcon.github.io/entropart/reference/BetaEntropy.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Estimate Shannon beta diversity
summary(BetaDiversity(Paracou618.MC, 1))
#> Neutral beta diversity of order 1 of metaCommunity Paracou618.MC 
#> with correction: Best 
#> 
#> Average diversity of the communities: 
#> [1] 1.436039
# Compare without correction
summary(BetaDiversity(Paracou618.MC, 1, Correction = "None"))
#> Neutral beta diversity of order 1 of metaCommunity Paracou618.MC 
#> with correction: None 
#> 
#> Average diversity of the communities: 
#> [1] 1.563888
# Estimate phylogenetic Shannon beta diversity
summary(BetaDiversity(Paracou618.MC, 1, Tree = Paracou618.Taxonomy) -> e)
#> HCDT beta diversity of order 1 of metaCommunity Paracou618.MC with correction: 
#> Best 
#> 
#> Phylogenetic or functional diversity was calculated according to the tree 
#> Paracou618.Taxonomy 
#> 
#> Diversity is normalized 
#> 
#> Average diversity of the communities: 
#> [1] 1.289993
```
