# Finds the original value (name of expression) of the argument of a function even in the case of embedded calls.

This function is used internally to find the name of arguments passed to
entropart functions such as
[`PhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)
that store them in their results.

## Usage

``` r
ArgumentOriginalName(x)
```

## Arguments

- x:

  Any argument.

## Value

The name of the argument.

## Details

The function searches the name of the argument in the parent frame of
the function and stops at the top.

## Author

BrodieG, \<https://stackoverflow.com/users/2725969/brodieg\>

Eric Marcon \<Eric.Marcon@agroparistech.fr\>
