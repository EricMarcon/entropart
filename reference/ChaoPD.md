# Phylogenetic Diversity of a Community

Calculates the phylogenetic diversity of order \\q\\ of a probability
vector.

## Usage

``` r
ChaoPD(Ps, q = 1, PhyloTree, Normalize = TRUE, Prune = FALSE, CheckArguments = TRUE)
```

## Arguments

- Ps:

  A probability vector, summing to 1.

- q:

  A number: the order of diversity. Default is 1.

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

The phylogenetic diversity is calculated following Chao *et al.* (2010).
The result is identical to the total diversity calculated by
[`PhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)
but it is much faster. A single value is returned instead of a
[`PhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)
object, and no bias correction is available.

The `Normalize` arguments allows calculating either \\^{q}\bar{D}(T)\\
(if `TRUE`) or \\^{q}PD(T)\\ if `FALSE`.

Diversity can be calculated for non ultrametric trees following Chao *et
al.* (2010) even though the meaning of the result is not so clear
(Leinster and Cobold, 2011).

## Value

A named number equal the diversity of the community. The name is "None"
to recall that no bias correction is available.

## References

Chao, A., Chiu, C.-H. and Jost, L. (2010). Phylogenetic diversity
measures based on Hill numbers. *Philosophical Transactions of the Royal
Society B* 365(1558): 3599-609.

Leinster, T. and Cobbold, C. (2011). Measuring diversity: the importance
of species similarity. *Ecology* 93(3): 477-489.

## See also

[`PhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md),
[`AllenH`](https://ericmarcon.github.io/entropart/reference/AllenH.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- Paracou618.MC$Ps

# Calculate the phylogenetic Simpson diversity of the plot
(ChaoPD(Paracou618.MC$Ps, 2, Paracou618.Taxonomy, Normalize=TRUE))
#>     None 
#> 24.61698 

# Calculate it using PhyloDiversity
# (more powerful but much slower if the tree has many periods)
PhyloDiversity(Paracou618.MC$Ps, 2, Paracou618.Taxonomy, Normalize=TRUE) -> phyD
summary(phyD)
#> alpha or gamma phylogenetic or functional diversity of order 2 of distribution 
#> -
#>  with correction: None
#> Phylogenetic or functional diversity was calculated according to the tree 
#> Paracou618.Taxonomy 
#> 
#> Diversity is normalized
#> 
#> Diversity equals: 24.61698 
```
