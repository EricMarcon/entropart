# Phylogenetic Diversity / Functional Diversity of a Community

Calculates Faith's PD / Petchey and Gaston' FD of a community described
by a probability vector and a phylogenetic / functional tree.

## Usage

``` r
PDFD(Ps, Tree, CheckArguments = TRUE)
```

## Arguments

- Ps:

  A probability vector, summing to 1.

- Tree:

  An object of class [`hclust`](https://rdrr.io/r/stats/hclust.html),
  "phylo" (see
  [`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)),
  [`phylog`](https://adeverse.github.io/ade4/reference/phylog.html) or
  [`PPtree`](https://ericmarcon.github.io/entropart/reference/PPtree.md).
  The tree must be ultrametric.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

PD and FD are defined as the total legth of the branches of the tree.

The probability vector is used to select branches: branches with
probability 0 are eliminated.

Bias correction requires the number of individuals to estimate sample
[`Coverage`](https://ericmarcon.github.io/entropart/reference/Coverage.md).

Use
[`bcPhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)`(Ps, 0, Tree)`
and choose the `Correction`.

## Value

A named number equal to the calculated diversity. The name is that of
the bias correction used.

## References

Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity.
*Biological Conservation* 61(1): 1-10.

Petchey, O. L. and Gaston, K. J. (2002). Functional diversity (FD),
species richness and community composition. *Ecology Letters* 5:
402-411.

## See also

[`bcPhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- Paracou618.MC$Ps
# Calculate the phylogenetic Shannon diversity of the plot
PDFD(Ps, Paracou618.Taxonomy)
#> None 
#>  395 
```
