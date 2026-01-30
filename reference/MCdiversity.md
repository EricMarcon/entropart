# Meta-Community diversity class.

Methods for objects of type "`MCdiversity`".

## Usage

``` r
is.MCdiversity(x)
# S3 method for class 'MCdiversity'
plot(x, ...)
# S3 method for class 'MCdiversity'
autoplot(object, col = "grey35",
          border = NA, ...)
# S3 method for class 'MCdiversity'
summary(object, ...)
```

## Arguments

- x:

  An object to be tested or plotted.

- object:

  A `MCdiversity` object to be summarized or plotted.

- col:

  The color used to fill the bars. See "Color Specification" in
  [`par`](https://rdrr.io/r/graphics/par.html).

- border:

  The color of the borders around the bars. See
  [`hist`](https://rdrr.io/r/graphics/hist.html).

- ...:

  Additional arguments to be passed to the generic methods.

## Value

Meta-community diversity objects are lists containing:

- MetaCommunity:

  The name of the
  [`MetaCommunity`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md)
  object containing inventory data.

- Type:

  The type of diversity ("alpha", "beta" or "gamma").

- Order:

  The order of diversity \\q\\.

- Correction:

  The estimation bias correction used to calculate diversity.

- Tree:

  The phylogenetic or functional tree used to calculate phylodiversity.

- Normalized:

  Logical. Indicates whether phylodiversity is normalized or
  proportional to the height of the tree.

- Weights:

  A vector containing the weights of communities.

- Communities:

  A vector containing the diversity of communities.

- Total:

  The total diversity.

`is.MCdiversity` returns `TRUE` if the object is of class `MCdiversity`.

`summary.MCdiversity` returns a summary of the object's value.
