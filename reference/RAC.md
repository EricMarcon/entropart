# Fit Distributions to Well-Known Rank Abundance Curves.

Observed distributions are fitted to classical RAC's.

## Usage

``` r
RAClnorm(Ns, CheckArguments = TRUE)
RACgeom(Ns, CheckArguments = TRUE)
RAClseries(Ns, CheckArguments = TRUE)
RACbstick(Ns, CheckArguments = TRUE)
```

## Arguments

- Ns:

  A numeric vector containing species abundances.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

`SpeciesDistribution` or integer vectors can be used to fit classical
rank-abundance curves (RAC) of classical distributions: `"RAClnorm"` for
log-normal (Preston, 1948), `"RAClseries"` for log-series (Fisher *et
al.*, 1943), `"RACgeom"` for geometric (Motomura, 1932) or `"RACbstick"`
for broken stick (MacArthur, 1957). method returns the estimated
parameters of the fitted distribution. The broken stick has no
parameter, so the maximum abundance is returned.

## Note

Fisher's alpha is estimated to fit the log-series distribution. The
estimation is done by the
[`fisher.alpha`](https://vegandevs.github.io/vegan/reference/diversity.html)
function of package `vegan`. It may differ substantially from the
estimation returned by `optimal.theta` from package `untb`.

## Value

A list (the parameters of distributions are returned only if the
distribution has been fit):

- Rank:

  A numeric vector. The ranks of species in the fitted RAC.

- Abundance:

  The abundance of species in the fitted RAC.

- mu:

  The expectation of the log-normal distribution

- sigma:

  The standard deviation of the log-normal distribution

- alpha:

  Fisher's alpha in the log-series distribution

- prob:

  The proportion of ressources taken by successive species in the
  geometric distribution

- max:

  The maximum abundance in the broken-stick distribution

## Author

Eric Marcon \<Eric.Marcon@agroparistech.fr\>, Bruno Herault
\<Bruno.Herault@cirad.fr\>

## References

Fisher R.A., Corbet A.S., Williams C.B. (1943) The Relation Between the
Number of Species and the Number of Individuals in a Random Sample of an
Animal Population. *Journal of Animal Ecology* 12: 42-58.

MacArthur R.H. (1957) On the Relative Abundance of Bird Species. *PNAS*
43(3): 293-295.

Motomura I. (1932) On the statistical treatment of communities.
*Zoological Magazine* 44: 379-383.

Preston, F.W. (1948). The commonness, and rarity, of species. *Ecology*
29(3): 254-283.

## See also

[`rgeom`](https://rdrr.io/r/stats/Geometric.html),
[`rlnorm`](https://rdrr.io/r/stats/Lognormal.html),
[`rCommunity`](https://ericmarcon.github.io/entropart/reference/rCommunity.md),
[`plot.SpeciesDistribution`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Fitted parameters
RACln <- RAClnorm(Ns)
RACln$mu
#> [1] 0.9519812
RACln$sigma
#> [1] 1.011383
RACgeom(Ns)$prob
#> [1] 0.01417454
RAClseries(Ns)$alpha
#> [1] 86.94339
RACbstick(Ns)$max
#> [1] 29.51418
```
