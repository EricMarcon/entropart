# Phylogenetic entropy or diversity.

Entropy or diversity against the height of the phylogenetic or
functional tree.

## Usage

``` r
is.PhyloValue(x)
# S3 method for class 'PhyloValue'
autoplot(object, xlab = expression(italic("T")), ylab = NULL, main = NULL,
          col = "black",
          lty = 1,
          lwd = 0.5,
          ...)
# S3 method for class 'PhyloValue'
plot(x, xlab = expression(italic("T")), ylab = NULL, main = NULL, ...)
# S3 method for class 'PhyloValue'
summary(object, ...)
```

## Arguments

- x:

  An object of class `PhyloValue`, including
  [`PhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)
  and
  [`PhyloEntropy`](https://ericmarcon.github.io/entropart/reference/PhyloEntropy.md)
  objects.

- xlab:

  The X axis label, "T" by default for Time.

- ylab:

  The Y axis label. if `NULL` (by default), "Entropy" or "Diversity" or
  nothing is chosen according to the object class.

- main:

  The main title of the plot. if `NULL` (by default), a default value is
  used.

- object:

  A `PhyloValue` object to be summarized.

- col:

  The color of the geom objects. See "Color Specification" in
  [`par`](https://rdrr.io/r/graphics/par.html).

- lty:

  The type of the lines. See
  [`lines`](https://rdrr.io/r/graphics/lines.html).

- lwd:

  The width of the lines. See
  [`lines`](https://rdrr.io/r/graphics/lines.html).

- ...:

  Additional arguments to be passed to
  [`plot`](https://rdrr.io/r/base/plot.html).

## Details

`PhyloValue` objects are the result of
[`PhyloApply`](https://ericmarcon.github.io/entropart/reference/PhyloApply.md).

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
#      and their taxonomy)
data(Paracou618)
# Calculate richness along the tree
# (Cuts are 1=species, 2=genus, 3=family)
summary(r <- PhyloApply(Paracou618.Taxonomy, FUN=Richness,
  NorP=Paracou618.MC$Ns, Normalize=TRUE))
#> Richness applied to Paracou618.MC$Ns along the tree: Paracou618.Taxonomy
#> 
#> Results are normalized
#> 
#> The average value is: 191
#> 
#> Values along the tree are:
#>   1   2   3 
#> 359 156  58 
autoplot(r)
```
