# Functional distances between pairs of species of Paracou field station plots 6 and 18, two 1-ha plots inventoried by the Bridge project.

This dataset is from Paracou field station, French Guiana, managed by
[Cirad](https://www.cirad.fr). Traits are detailed in Marcon and Herault
(2014), the distance matrix was built following Paine et al. (2011).

## Usage

``` r
data(Paracou618)
```

## Format

An object of class [`dist`](https://rdrr.io/r/stats/dist.html).

## Source

Permanent data census of Paracou.

## References

Gourlet-Fleury, S., Guehl, J. M. and Laroussinie, O., Eds. (2004).
*Ecology & management of a neotropical rainforest. Lessons drawn from
Paracou, a long-term experimental research site in French Guiana*.
Paris, Elsevier.

Marcon, E., Herault, B. (2015). Decomposing Phylodiversity. *Methods in
Ecology and Evolution* 6(3): 333-339.

Paine, C. E. T., Baraloto, C., Chave, J., and Herault, B. (2011).
Functional traits of individual trees reveal ecological constraints on
community assembly in tropical rain forests. *Oikos*, 120(5), 720-727.

## Examples

``` r
data(Paracou618)
plot(density(Paracou618.dist, from=0), main="Distances between species")
```
