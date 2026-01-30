# Checks the arguments of a function of the package entropart

This function is used internally to verify that arguments passed to
entropart functions such as
[`PhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)
are correct.

## Usage

``` r
CheckentropartArguments()
```

## Value

Returns `TRUE` or stops if a problem is detected.

## Details

The function compares the arguments passed to its parent function to the
type they should be and performs some extra tests (*e.g.* probabilities
must be positive and sum to 1). It stops if an argument is not correct.
