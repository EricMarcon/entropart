# Exponential of order q

Calculates the deformed exponential of order \\q\\.

## Usage

``` r
expq(x, q)
expq.CommunityProfile(Profile)
```

## Arguments

- x:

  A numeric vector.

- Profile:

  A
  [`CommunityProfile`](https://ericmarcon.github.io/entropart/reference/CommunityProfile.md).

- q:

  A number.

## Details

The deformed exponential is defined as \\(x(1-q)+1)^{\frac{1}{(1-q)}}\\.

For \\q\>1\\, \\\ln_q{(+\infty)}=\frac{1}{(q-1)}\\ so \\\exp_q{(x)}\\ is
not defined for \\x\>\frac{1}{(q-1)}\\.

`expq.CommunityProfile` calculates the deformed exponential of a
[`CommunityProfile`](https://ericmarcon.github.io/entropart/reference/CommunityProfile.md).
Its `$x` item (the order of diversity) is kept unchanged whilst other
items are set to their exponential of order `$x`. Thus, an entropy
profile is transformed into a diversity profile.

## Value

A vector of the same length as `x` containing the transformed values or
a
[`CommunityProfile`](https://ericmarcon.github.io/entropart/reference/CommunityProfile.md).

## References

Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014).
Generalization of the partitioning of Shannon diversity. *PLOS One*
9(3): e90289.

Tsallis, C. (1994). What are the numbers that experiments provide?
*Quimica Nova* 17(6): 468-471.

## See also

`expq`

## Examples

``` r
curve(exp(x), -5, 0, lty=3)  
curve(expq(x, 2), -5, 0, lty=2, add=TRUE)
curve(expq(x, 3), -5, 0, lty=1, add=TRUE)
legend("topleft", legend = c("exp(x)", "exp2(x)", "exp3(x)"), lty = c(1, 2, 3), inset=0.02)
```
