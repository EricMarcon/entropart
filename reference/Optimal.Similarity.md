# Optimal scale parameter to transform a distance matrix into a similarity matrix

Calculates the scale parameter \\u\\ that maximizes the variance of the
similarity matrix \\exp(-u\*DistanceMatrix)\\.

## Usage

``` r
Optimal.Similarity(Distance, CheckArguments = TRUE)
```

## Arguments

- Distance:

  A distance matrix, *i.e.* a square matrix with zeros on its diagonal
  or a [`dist`](https://rdrr.io/r/stats/dist.html) object.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

The similarity matrix used by
[`Dqz`](https://ericmarcon.github.io/entropart/reference/Dqz.md)) can be
optimized following Marcon *et al.* (2014) such that the variance of
similarities between pairs of species is maximized. See the example.

## Value

A list:

- u:

  The optimal scale \\u\\.

- Matrix:

  The optimal similarity matrix `Z`.

## References

Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of
similarity-based diversity and its bias correction. *HAL*
hal-00989454(version 3).

## See also

[`Dqz`](https://ericmarcon.github.io/entropart/reference/Dqz.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Prepare the similarity matrix. The square root of Paracou618.dist is euclidean.
optimal <- Optimal.Similarity(sqrt(Paracou618.dist))
# Optimal scale
optimal$u
#> [1] 2.348622
# Calculate diversity of order 2
bcDqz(Paracou618.MC$Ns, 2, optimal$Matrix)
#>     Best 
#> 1.558321 
```
