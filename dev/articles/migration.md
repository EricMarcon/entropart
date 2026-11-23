# Migration

*entropart* is being superseded by the
[*divent*](https://ericmarcon.github.io/divent/) package. It will be
maintained, i.e. bugs will be corrected and updates required by CRAN
will be made, but no new development is expected.

## General changes

Object names follow the Pascal-case style in *entropart*, e.g.
[`DivPart()`](https://ericmarcon.github.io/entropart/dev/reference/DivPart.md).
*divent* follows the recommendations of the tidyverse, e.g. snake-case
names with prefixes,
e.g. [`div_part()`](https://ericmarcon.github.io/divent/reference/div_part.html).

Functions’ values are often lists in *entropart*, vs. tibbles in
*divent*, where argument `as_numeric = TRUE` allows obtaining numeric
outputs from many functions.

Abundance or probability data are named vectors in *entropart*, declared
as `AbdVector` or `ProbaVector` objects by the
[`as.AbdVector()`](https://ericmarcon.github.io/entropart/dev/reference/SpeciesDistribution.md)
or
[`as.ProbaVector()`](https://ericmarcon.github.io/entropart/dev/reference/SpeciesDistribution.md)function.
In *divent*, numeric vectors are used directly or converted to tibbles
of class `abundances` or `probabilities` by functions
[`as_abundances()`](https://ericmarcon.github.io/divent/reference/species_distribution.html)
or
[`as_probabilities()`](https://ericmarcon.github.io/divent/reference/species_distribution.html).

The following table shows the correspondence between *entropart* and
*divent* functions.

## Function names

[TABLE]

## Example

``` r

library("entropart")
library("divent")
# Random community, lognormal, 10 individuals, 3 species 
entropart::rCommunity(
  n = 1,
  size = 10, 
  S = 3, 
  Distribution = "lnorm"
)
```

    ## [1] 3 4 3
    ## attr(,"class")
    ## [1] "AbdVector"           "SpeciesDistribution" "integer"

``` r

# Equivalent to
divent::rcommunity(
  n = 1, 
  size = 10, 
  species_number = 3, 
  distribution = "lnorm"
)
```

    ## # A tibble: 1 × 5
    ##   site    weight  sp_1  sp_2  sp_3
    ##   <chr>    <dbl> <int> <int> <int>
    ## 1 lnorm_1     10     2     7     1

``` r

# use as.numeric() to produce a vector rather than a tibble
divent::rcommunity(
  n = 1, 
  size = 10, 
  species_number = 3, 
  distribution = "lnorm"
) |>
  as.numeric()
```

    ## sp_1 sp_2 sp_3 
    ##    8    2    0
