# Grassberger's expectation of n^q

Expected value of \\n^q\\ when \\n\\ follows a Poisson law.

## Usage

``` r
Enq(n, q)
```

## Arguments

- n:

  A positive integer vector.

- q:

  A positive number.

## Details

The expectation of \\n^q\\ when \\n\\ follows a Poisson ditribution has
been derived by Grassberger (1988).

## Value

A vector of the same length as `n` containing the transformed values.

## References

Grassberger, P. (1988). Finite sample corrections to entropy and
dimension estimates. Physics Letters A 128(6-7): 369-373.

## Note

The function is computed using the
[`beta`](https://rdrr.io/r/base/Special.html).function.

Its value is 0 for \\n-q+1\<0\\.

## Examples

``` r
# Compare
n <- c(2,3)
Enq(n, q=2)
#> [1] 2 6
# with
n^2
#> [1] 4 9

# Result is 1
Enq(n, q=0)
#> [1] 1 1
# Result is 0
Enq(n, q=5)
#> [1] -7.796344e-15  5.847258e-15
```
