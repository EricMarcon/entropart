# Generalized Kullback-Leibler divergence

Calculates the generalized Kullback-Leibler divergence between an
observed and an expected probability distribution.

## Usage

``` r
KLq(Ps, Pexp, q = 1, CheckArguments = TRUE)
```

## Arguments

- Ps:

  The observed probability vector.

- Pexp:

  The expected probability vector.

- q:

  A number: the order of entropy. Default is 1.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

The generalized Kullback-Leibler divergence (Borland et al., 1998)
converges to the Kullback-Leibler divergence (Kullback and Leibler,
1951) when \\q\\ tends to 1. It is used to calculate the generalized
beta entropy (Marcon et al., 2014).

## Value

A number equal to the generalized Kullback-Leibler divergence between
the probability distributions.

## References

Borland, L., Plastino, A. R. and Tsallis, C. (1998). Information gain
within nonextensive thermostatistics. *Journal of Mathematical Physics*
39(12): 6490-6501.

Kullback, S. and Leibler, R. A. (1951). On Information and Sufficiency.
The Annals of Mathematical *Statistics* 22(1): 79-86.

Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014).
Generalization of the partitioning of Shannon diversity. *PLOS One*
9(3): e90289.

## See also

[`TsallisBeta`](https://ericmarcon.github.io/entropart/reference/TsallisBeta.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ps is the vector of probabilities
Ps <- Paracou618.MC$Ps
# Probability distribution of the first plot
Ps1 <- Paracou618.MC$Psi[, 1]
# Divergence of order 2 between the first plot and the whole forest
KLq(Ps1, Ps, 2)
#> [1] 0.01004705
```
