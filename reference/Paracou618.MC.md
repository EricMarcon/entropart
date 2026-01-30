# Paracou field station plots 6 and 18, two 1-ha plots inventoried by the Bridge project.

This dataset is from Paracou field station, French Guiana, managed by
[Cirad](https://www.cirad.fr).

## Usage

``` r
data(Paracou618)
```

## Format

An object of class
[`MetaCommunity`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md)
made of two communities and 425 species.

## Source

Permanent data census of Paracou and Marcon et al. (2012).

## References

Gourlet-Fleury, S., Guehl, J. M. and Laroussinie, O., Eds. (2004).
*Ecology & management of a neotropical rainforest. Lessons drawn from
Paracou, a long-term experimental research site in French Guiana*.
Paris, Elsevier.

Marcon, E., F. Puech, et al. (2012). Characterizing the relative spatial
structure of point patterns. *International Journal of Ecology*
2012(Article ID 619281): 11.

## Examples

``` r
data(Paracou618)
summary(Paracou618.MC)
#> Meta-community (class 'MetaCommunity') made of 1124 individuals in 2 
#> communities and 425 species. 
#> 
#> Its sample coverage is 0.92266748426447 
#> 
#> Community weights are: 
#> [1] 0.5720641 0.4279359
#> Community sample numbers of individuals are: 
#> P006 P018 
#>  643  481 
#> Community sample coverages are: 
#>      P006      P018 
#> 0.8943859 0.8463782 
```
