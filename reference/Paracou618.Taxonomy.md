# Taxonomy (Family - Genus - Species) of Paracou field station plots 6 and 18, two 1-ha plots inventoried by the Bridge project.

This dataset is from Paracou field station, French Guiana, managed by
[Cirad](https://www.cirad.fr).

## Usage

``` r
data(Paracou618)
```

## Format

An object of class "phylo" (see
[`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)) containing a
taxonomy.

## Source

Permanent data census of Paracou.

## References

Gourlet-Fleury, S., Guehl, J. M. and Laroussinie, O., Eds. (2004).
*Ecology & management of a neotropical rainforest. Lessons drawn from
Paracou, a long-term experimental research site in French Guiana*.
Paris, Elsevier.

## Examples

``` r
data(Paracou618)
plot(Paracou618.Taxonomy, type="fan", show.tip.label=FALSE)
```
