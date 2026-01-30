# Similarity-based entropy of a community

Calculates the entropy of order \\q\\ of a probability vector according
to a similarity matrix.

## Usage

``` r
Hqz(NorP, q = 1, Z = diag(length(NorP)),  ...)
bcHqz(Ns, q = 1, Z = diag(length(Ns)), Correction = "Best", SampleCoverage = NULL, 
  CheckArguments = TRUE)
# S3 method for class 'ProbaVector'
Hqz(NorP, q = 1, Z = diag(length(NorP)), 
  ..., CheckArguments = TRUE, Ps = NULL)
# S3 method for class 'AbdVector'
Hqz(NorP, q = 1, Z = diag(length(NorP)), Correction = "Best",
  ..., CheckArguments = TRUE, Ns = NULL)
# S3 method for class 'integer'
Hqz(NorP, q = 1, Z = diag(length(NorP)), Correction = "Best", 
  ..., CheckArguments = TRUE, Ns = NULL)
# S3 method for class 'numeric'
Hqz(NorP, q = 1, Z = diag(length(NorP)), Correction = "Best", 
  ..., CheckArguments = TRUE, Ps = NULL, Ns = NULL)
```

## Arguments

- Ps:

  A probability vector, summing to 1.

- Ns:

  A numeric vector containing species abundances.

- NorP:

  A numeric vector, an integer vector, an abundance vector
  ([`AbdVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md))
  or a probability vector
  ([`ProbaVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md)).
  Contains either abundances or probabilities.

- q:

  A number: the order of entropy. Default is 1.

- Z:

  A relatedness matrix, *i.e.* a square matrix whose terms are all
  positive, strictly positive on the diagonal. Generally, the matrix is
  a similarity matrix, *i.e.* the diagonal terms equal 1 and other terms
  are between 0 and 1. Default is the identity matrix to calculate
  neutral entropy.

- Correction:

  A string containing one of the possible corrections: `"None"` (no
  correction), `"ChaoShen"`, `"MarconZhang"` or `"Best"`, the default
  value. The `"MarconZhang"` correction assumes a similarity matrix.

- SampleCoverage:

  The sample coverage of `Ns` calculated elsewhere. Used to calculate
  the gamma diversity of meta-communities, see details.

- ...:

  Additional arguments. Unused.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

Entropy is calculated following Leinster and Cobbold (2012) after
Ricotta and Szeidl (2006): it is the entropy of order `q` of the
community, using species ordinariness as the information function.

A similarity matrix is used (as for
[`Dqz`](https://ericmarcon.github.io/entropart/reference/Dqz.md)), not a
distance matrix as in Ricotta and Szeidl (2006). See the example.

Bias correction requires the number of individuals. Use `bcHqz` and
choose the `Correction`. Correction techniques are from Marcon *et al.*
(2014).

Currently, the `"Best"` correction is the max value of `"ChaoShen"` and
`"MarconZhang"`.

The functions are designed to be used as simply as possible. `Hqz` is a
generic method. If its first argument is an abundance vector, an integer
vector or a numeric vector which does not sum to 1, the bias corrected
function `bcHqz` is called. Explicit calls to `bcHqz` (with bias
correction) or to `Hqz.ProbaVector` (without correction) are possible to
avoid ambiguity. The `.integer` and `.numeric` methods accept `Ps` or
`Ns` arguments instead of `NorP` for backward compatibility.

The size of a metacommunity (see
[`MetaCommunity`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md))
is unknown so it has to be set according to a rule which does not ensure
that its abundances are integer values. Then, classical bias-correction
methods do not apply. Providing the `SampleCoverage` argument allows
applying the `"ChaoShen"` correction to estimate quite well the entropy.
[`DivPart`](https://ericmarcon.github.io/entropart/reference/DivPart.md)
and
[`GammaEntropy`](https://ericmarcon.github.io/entropart/reference/GammaEntropy.md)
functions use this tweak.

## Value

A named number equal to the calculated entropy. The name is that of the
bias correction used.

## References

Leinster, T. and Cobbold, C. (2012). Measuring diversity: the importance
of species similarity. *Ecology* 93(3): 477-489.

Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of
similarity-based diversity and its bias correction. *HAL*
hal-00989454(version 3).

Ricotta, C. and Szeidl, L. (2006). Towards a unifying approach to
diversity measures: Bridging the gap between the Shannon entropy and
Rao's quadratic index. *Theoretical Population Biology* 70(3): 237-243.

## See also

[`Dqz`](https://ericmarcon.github.io/entropart/reference/Dqz.md),
[`PhyloEntropy`](https://ericmarcon.github.io/entropart/reference/PhyloEntropy.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Prepare the similarity matrix
DistanceMatrix <- as.matrix(EightSpTree$Wdist^2/2)
# Similarity can be 1 minus normalized distances between species
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
# Calculate diversity of order 2
Ps <- EightSpAbundance/sum(EightSpAbundance)
Hqz(Ps, 2, Z)
#>      None 
#> 0.6033223 
# Equal to normalized Rao quadratic entropy when q=2
Rao(Ps, EightSpTree)/max(DistanceMatrix)
#>      None 
#> 0.6033223 
# But different from PhyloEntropy for all other q, e.g. 1
Hqz(Ps, 1, Z)
#>      None 
#> 0.9348275 
summary(PhyloEntropy(Ps, 1, EightSpTree))
#> alpha or gamma phylogenetic or functional entropy of order 1 of distribution -
#>  with correction: None
#> Phylogenetic or functional entropy was calculated according to the tree 
#> EightSpTree 
#> 
#> Entropy is normalized
#> 
#> Entropy equals: 1.034531 
```
