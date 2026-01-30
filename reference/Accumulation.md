# Diversity accumulation.

Diversity and Entropy Accumulation Curves represent the accumulation of
entropy with respect to the sample size.

## Usage

``` r
as.AccumCurve(x, y, low = NULL, high = NULL)
is.AccumCurve(x)
EntAC(Ns, q = 0, n.seq = seq_len(sum(Ns)), PCorrection="Chao2015", Unveiling="geom",
  RCorrection="Rarefy", NumberOfSimulations = 0, Alpha = 0.05, 
  ShowProgressBar = TRUE, CheckArguments = TRUE)
DivAC(Ns, q = 0, n.seq = seq_len(sum(Ns)), PCorrection="Chao2015", Unveiling="geom",
  RCorrection="Rarefy", NumberOfSimulations = 0, Alpha = 0.05, 
  ShowProgressBar = TRUE, CheckArguments = TRUE)
# S3 method for class 'AccumCurve'
plot(x, ..., main = NULL, 
         xlab = "Sample Size", ylab = NULL, ylim = NULL,
         LineWidth = 2, ShadeColor = "grey75", BorderColor = "red")
# S3 method for class 'AccumCurve'
autoplot(object, ..., main = NULL, 
         xlab = "Sample Size", ylab = NULL, 
         ShadeColor = "grey75", alpha = 0.3, BorderColor = "red",
         col = "black",
         lty = 1,
         lwd = 0.5)
```

## Arguments

- x:

  An object. A numeric vector in `as.AccumCurve`.

- object:

  An object.

- y:

  A numeric vector.

- low:

  A numeric vector.

- high:

  A numeric vector.

- Ns:

  A numeric vector containing species abundances.

- q:

  A number: the order of diversity. Default is 1.

- n.seq:

  A sequence of numbers. Accumulation will be calculated at each value.

- PCorrection:

  A string containing one of the possible corrections to estimate a
  probability distribution in
  [`as.ProbaVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md):
  `"Chao2015"` is the default value. Used only for extrapolation and `q`
  different from 0, 1, 2.

- Unveiling:

  A string containing one of the possible unveiling methods to estimate
  the probabilities of the unobserved species in
  [`as.ProbaVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md):
  `"geom"` (geometric: the unobserved species distribution is geometric)
  is the default value. Used only for extrapolation and `q` different
  from 0, 1, 2.

- RCorrection:

  A string containing a correction recognized by
  [`Richness`](https://ericmarcon.github.io/entropart/reference/Richness.md)
  to evaluate the total number of species in
  [`as.ProbaVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md).
  `"Rarefy"` is the default value to estimate the number of species such
  that the entropy of the asymptotic distribution rarefied to the
  observed sample size equals the observed entropy of the data. Used
  only for extrapolation and `q` different from 0, 1, 2. If `q` is 0
  (extrapolation of richness), `"Rarefy"` is taken for `"Jackknife"`.

- NumberOfSimulations:

  The number of Simulations to build confidence intervals.

- Alpha:

  The risk level, 5% by default.

- ...:

  Additional arguments to be passed to
  [`plot`](https://rdrr.io/r/base/plot.html). Unused elsewhere.

- main:

  The main title of the plot. if `NULL` (by default), there is no title.

- xlab:

  The X axis label, "Rank" by default.

- ylab:

  The Y axis label. if `NULL` (by default), "Probability" or "Abundance"
  is chosen according to the object class.

- ylim:

  The interval of y values plotted.

- LineWidth:

  The width of the line that represents the actual profile.

- ShadeColor:

  The color of the shaded confidence envelope.

- BorderColor:

  The color of the bounds of the confidence envelope.

- alpha:

  Opacity of the confidence enveloppe, between 0 and 1.

- col:

  The color of the geom objects. See "Color Specification" in
  [`par`](https://rdrr.io/r/graphics/par.html).

- lty:

  The type of the lines. See
  [`lines`](https://rdrr.io/r/graphics/lines.html).

- lwd:

  The width of the lines. See
  [`lines`](https://rdrr.io/r/graphics/lines.html).

- ShowProgressBar:

  If `TRUE` (default), a progress bar is shown.

- CheckArguments:

  Logical; if `TRUE`, the function arguments are verified. Should be set
  to `FALSE` to save time when the arguments have been checked
  elsewhere.

## Details

`DivAC` or `EntAC` estimate the diversity or entropy accumulation curve
of a distribution. See
[`Tsallis`](https://ericmarcon.github.io/entropart/reference/Tsallis.md)
for details about the computation of entropy at each level of
interpolation and extrapolation. In accumulation curves, extrapolation
if done by estimating the asymptotic distribution of the community and
estimating entropy at different levels by interpolation. The asymptotic
richess is adjusted so that the extrapolated part of the accumulation
joins the observed value at the sample size.

`AccumCurve` objects include `EntAC` and `DivAC` objects for entropy and
diversity accumulation. They generalize the classical Species
Accumulation Curves (SAC) which are diversity accumulation of order
\$q=0\$.

`as.AccumCurve` transforms two vectors (where `x` is the sammple size
and `y` the accumulation) into an object of class `AccumCurve`.

`AccumCurve` objects can be plotted with either `plot` or
`autoplot methods`.

## Value

A `DivAC` or an `EntAC` object. Both are `AccumCurve` objects, which are
a list:

- x:

  The sample size.

- y:

  The value of entropy or diversity.

- low:

  The lower bound of the confidence envelope of the estimation.

- high:

  The upper bound of the confidence envelope of the estimation.

Attibutes `"Size"` and `"Value"` contain the actual sample size and the
corresponding diversity or entropy.

`AccumCurve` objects can be summarized and plotted.

## References

Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H.,
Colwell, R. K., Ellison, A. M (2014). Rarefaction and extrapolation with
Hill numbers: A framework for sampling and estimation in species
diversity studies. *Ecological Monographs*, 84(1): 45-67.

## See also

[`Tsallis`](https://ericmarcon.github.io/entropart/reference/Tsallis.md),
[`Diversity`](https://ericmarcon.github.io/entropart/reference/Diversity.md)

## Examples

``` r
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Accumulation curve of Simpson's diversity
autoplot(DivAC(Ns, q=2))
#> 
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the entropart package.
#>   Please report the issue at <https://github.com/EricMarcon/entropart/issues/>.
```
