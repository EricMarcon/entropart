# Phylogenetic Entropy of a Community

Calculates the phylogenetic diversity of order \\q\\ of a probability
vector.

## Usage

``` r
AllenH(Ps, q = 1, PhyloTree, Normalize = TRUE, Prune = FALSE, CheckArguments = TRUE)
```

## Arguments

- Ps:

  A probability vector, summing to 1.

- q:

  A number: the order of entropy. Default is 1.

- PhyloTree:

  An object of class [`hclust`](https://rdrr.io/r/stats/hclust.html),
  "phylo" (see
  [`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)),
  [`phylog`](https://adeverse.github.io/ade4/reference/phylog.html) or
  [`PPtree`](https://ericmarcon.github.io/entropart/reference/PPtree.md).
  The tree is not necessarily ultrametric.

- Normalize:

  If `TRUE` (default), diversity is not affected by the height of the
  tree.  
  If `FALSE`, it is proportional to the height of the tree.

- Prune:

  What to do when somes species are in the tree but not in `Ps`?  
  If `TRUE`, the tree is pruned to keep species of `Ps` only. The height
  of the tree may be changed if a pruned branch is related to the
  root.  
  If `FALSE` (default), species with probability 0 are added in `Ps`.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

The phylogenetic entropy is calculated following Allen *et al.* (2009)
for order \\q=1\\ and Leinster and Cobold (2011) for other orders.The
result is identical to the total entropy calculated by
[`PhyloEntropy`](https://ericmarcon.github.io/entropart/reference/PhyloEntropy.md)
but it is much faster. A single value is returned instead of a
[`PhyloEntropy`](https://ericmarcon.github.io/entropart/reference/PhyloEntropy.md)
object, and no bias correction is available.

The `Normalize` argument allows normalizing entropy by the height of the
tree, similarly to
[`ChaoPD`](https://ericmarcon.github.io/entropart/reference/ChaoPD.md).

Diversity can be calculated for non ultrametric trees following Leinster
and Cobold (2011) even though the meaning of the result is not so clear.

## Value

A named number equal the entropy of the community. The name is "None" to
recall that no bias correction is available.

## References

Allen, B., Kon, M. and Bar-Yam, Y. (2009). A New Phylogenetic Diversity
Measure Generalizing the Shannon Index and Its Application to
Phyllostomid Bats. *American Naturalist* 174(2): 236-243.

Leinster, T. and Cobbold, C. (2011). Measuring diversity: the importance
of species similarity. *Ecology* 93(3): 477-489.

## See also

[`PhyloEntropy`](https://ericmarcon.github.io/entropart/reference/PhyloEntropy.md),
[`ChaoPD`](https://ericmarcon.github.io/entropart/reference/ChaoPD.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ns)

# Calculate the phylogenetic Shannon diversity of the plot
AllenH(Ps, 1, Paracou618.Taxonomy, Normalize=TRUE)
#>     None 
#> 3.919941 

# Calculate it using PhyloEntropy: more powerful but much slower is the tree has many periods
PhyloEntropy(Ps, 1, Paracou618.Taxonomy, Normalize=TRUE) -> phyE
summary(phyE)
#> alpha or gamma phylogenetic or functional entropy of order 1 of distribution -
#>  with correction: None
#> Phylogenetic or functional entropy was calculated according to the tree 
#> Paracou618.Taxonomy 
#> 
#> Entropy is normalized
#> 
#> Entropy equals: 3.919941 
```
