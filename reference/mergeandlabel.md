# Merge community data

Merge two dataframes containing species abundances.

## Usage

``` r
mergeandlabel(x, y)
```

## Arguments

- x:

  A dataframe containing species abundances

- y:

  A dataframe containing species abundances

## Details

Row names must contain species. The new dataframe keeps all species with
abundance zero for species that were not found in original communities.

## Value

A dataframe containing all lines and columns of the merged dataframes.

## See also

[`MergeC`](https://ericmarcon.github.io/entropart/reference/MergeMC.md),
[`MergeMC`](https://ericmarcon.github.io/entropart/reference/MergeMC.md),
[`ShuffleMC`](https://ericmarcon.github.io/entropart/reference/MergeMC.md)
