# Preprocessing of a phylogenetic tree

Calculates cuts and intervals of a phylogenetic tree and make it
available both in [`hclust`](https://rdrr.io/r/stats/hclust.html) and
"phylo" (see [`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html))
formats.

## Usage

``` r
Preprocess.Tree(Tree)
```

## Arguments

- Tree:

  An object of class [`hclust`](https://rdrr.io/r/stats/hclust.html),
  "phylo" (see
  [`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)) or
  [`phylog`](https://adeverse.github.io/ade4/reference/phylog.html). The
  tree must be ultrametric.

## Details

This is an internal function. It is used by
[`PhyloApply`](https://ericmarcon.github.io/entropart/reference/PhyloApply.md)
to obtain values of intervals.

## Value

An object of class
[`PPtree`](https://ericmarcon.github.io/entropart/reference/PPtree.md).
