# Reduced-bias alpha diversity of a metacommunity

Calculates the eeduced-bias total alpha diversity of order \\q\\ of
communities.

## Usage

``` r
AlphaDiversity(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, 
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
[`AlphaEntropy`](https://ericmarcon.github.io/entropart/reference/AlphaEntropy.md)
and transformed into diversity.

## Value

An
[`MCdiversity`](https://ericmarcon.github.io/entropart/reference/MCdiversity.md)
object containing diversity values of each community and of the
metacommunity.

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

[`AlphaEntropy`](https://ericmarcon.github.io/entropart/reference/AlphaEntropy.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Calculate Simpson alpha diversity
summary(AlphaDiversity(Paracou618.MC, 2))
#> Neutral alpha diversity of order 2 of metaCommunity Paracou618.MC 
#> with correction: Best 
#> 
#> Diversity of communities: 
#>     P006     P018 
#> 38.78158 53.97842 
#> Average diversity of the communities: 
#> [1] 44.09399
# Compare without correction
summary(AlphaDiversity(Paracou618.MC, 2, Correction = "None"))
#> Neutral alpha diversity of order 2 of metaCommunity Paracou618.MC 
#> with correction: None 
#> 
#> Diversity of communities: 
#>     P006     P018 
#> 36.65653 48.67684 
#> Average diversity of the communities: 
#> [1] 40.98793
# Estimate phylogenetic Simpson alpha diversity
summary(AlphaDiversity(Paracou618.MC, 2, Tree = Paracou618.Taxonomy) -> e)
#> HCDT alpha diversity of order 2 of metaCommunity Paracou618.MC with correction: 
#> Best 
#> 
#> Phylogenetic or functional diversity was calculated according to the tree 
#> Paracou618.Taxonomy 
#> 
#> Diversity is normalized 
#> 
#> Diversity of communities: 
#>     P006     P018 
#> 17.16423 25.49632 
#> Average diversity of the communities: 
#> [1] 19.95488
plot(e)  
```
