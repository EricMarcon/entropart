# Meta-Community entropy class.

Methods for objects of type "`MCentropy`".

## Usage

``` r
is.MCentropy(x)
# S3 method for class 'MCentropy'
plot(x, ...)
# S3 method for class 'MCentropy'
autoplot(object, col = "grey35",
          border = NA, ...)
# S3 method for class 'MCentropy'
summary(object, ...)
```

## Arguments

- x:

  An object to be tested or plotted.

- object:

  A `MCentropy` object to be summarized or plotted.

- col:

  The color used to fill the bars. See "Color Specification" in
  [`par`](https://rdrr.io/r/graphics/par.html).

- border:

  The color of the borders around the bars. See
  [`hist`](https://rdrr.io/r/graphics/hist.html).

- ...:

  Additional arguments to be passed to the generic methods.

## Value

Meta-community entropy objects are lists containing:

- MetaCommunity:

  The name of the
  [`MetaCommunity`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md)
  object containing inventory data.

- Method:

  The method used to calculate entropy ("HCDT", "Similarity-based").

- Type:

  The type of entropy ("alpha", "beta" or "gamma").

- Order:

  The order of entropy \\q\\.

- Correction:

  The estimation bias correction used to calculate entropy.

- Tree:

  The phylogenetic or functional tree used to calculate phyloentropy.

- Normalized:

  Logical. Indicates whether phyloentropy is normalized or proportional
  to the height of the tree.

- Z:

  The matrix used to calculate similarity-based entropy.

- Weights:

  A vector containing the weights of communities.

- Communities:

  A vector containing the entropy of communities.

- Total:

  The total entropy.

`is.MCentropy` returns `TRUE` if the object is of class `MCentropy`.

`summary.MCentropy` returns a summary of the object's value.
