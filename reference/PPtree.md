# Preprocessed Trees.

Methods for objects of type "`PPtree`".

## Usage

``` r
is.PPtree(x)
# S3 method for class 'PPtree'
plot(x, ...)
```

## Arguments

- x:

  An object to be tested or plotted

- ...:

  Additional arguments to be passed to the generic methods.

## Value

An object of class **PPtree** is a list:

- phyTree:

  A "phylo" (see
  [`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)) tree

- hTree:

  A [`hclust`](https://rdrr.io/r/stats/hclust.html) tree

- Height:

  The height of the tree, that is to say the distance between root and
  leaves

- Cuts:

  A vector. Cut times of the tree (the distance from nodes to leaves)

- Intervals:

  A vector. The lengths of intervals between cuts

`is.PPtree` returns `TRUE` if the object is of class `PPtree`.

`plot.PPtree` plots it.

## Note

Versions up to 1.3 contained a
[`phylog`](https://adeverse.github.io/ade4/reference/phylog.html) tree,
now deprecated in `ade4`. A "phylo" (see
[`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)) tree is now
used.

See the dedicated vignette
(`vignette("Phylogenies", package="entropart")`) for more details.

## Examples

``` r
data(Paracou618)
# Preprocess a phylog object
ppt <- Preprocess.Tree(EightSpTree)
# Is it a preprocessed tree?
is.PPtree(ppt)
#> [1] TRUE
# Plot it
plot(ppt, hang=-1)

# Alternative plot
ade4::radial.phylog(EightSpTree)
```
