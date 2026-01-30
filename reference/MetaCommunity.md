# Metacommunity class

Methods for objects of type "`MetaCommunity`".

## Usage

``` r
MetaCommunity(Abundances, Weights = rep(1, ncol(Abundances)))
is.MetaCommunity(x)
# S3 method for class 'MetaCommunity'
summary(object, ...)
# S3 method for class 'MetaCommunity'
plot(x, ...)
```

## Arguments

- Abundances:

  A dataframe containing the number of observations (lines are species,
  columns are communities). The first column of the dataframe may
  contain the species names.

- Weights:

  A vector of positive numbers equal to community weights or a dataframe
  containing a vector named `Weights`. It does not have to be
  normalized. Weights are equal by default.

- x:

  An object to be tested or plotted.

- object:

  A `MetaCommunity` object to be summarized.

- ...:

  Additional arguments to be passed to the generic methods.

## Details

In the entropart package, individuals of different "species" are counted
in several "communities" which are agregated to define a
"metacommunity".

This is a naming convention, which may correspond to plots in a forest
inventory or any data organized the same way.

Alpha and beta entropies of communities are summed according to
`Weights` and the probability to find a species in the metacommunity is
the weighted average of probabilities in communities.

The simplest way to import data is to organize it into two text files.
The first file should contain abundance data: the first column named
`Species` for species names, and a column for each community.

The second file should contain the community weights in two columns. The
first one, named `Communities` should contain their names and the second
one, named `Weights`, their weights.

Files can be read and data imported by code such as:

    Abundances <- read.csv(file="Abundances.csv", row.names = 1)
    Weights <- read.csv(file="Weights.csv")
    MC <- MetaCommunity(Abundances, Weights)

## Value

An object of class **MetaCommunity** is a list:

- Nsi:

  A matrix containing abundance data, species in line, communities in
  column.

- Ns:

  A vector containing the number of individuals of each species.

- Ni:

  A vector containing the number of individuals of each community.

- N:

  The total number of individuals.

- Psi:

  A matrix whose columns are the probability vectors of communities
  (each of them sums to 1).

- Wi:

  A vector containing the normalized community weights (sum to 1).

- Ps:

  A vector containing the probability vector of the metacommunity.

- Nspecies:

  The number of species.

- Ncommunities:

  The number of communities.

- SampleCoverage:

  The sample coverage of the metacommunity.

- SampleCoverage.communities:

  A vector containing the sample coverages of each community.

`is.MetaCommunity` returns `TRUE` if the object is of class
`MetaCommunity`.

`summary.MetaCommunity` returns a summary of the object's value.

`plot.MetaCommunity` plots it.

## Examples

``` r
# Use BCI data from vegan package
if (require(vegan, quietly = TRUE)) {
  # Load BCI data (number of trees per species in each 1-ha plot of a tropical forest)
  data(BCI)
  # BCI dataframe must be transposed (its lines are plots, not species)
  BCI.df <- as.data.frame(t(BCI))
  # Create a metacommunity object from a matrix of abundances and a vector of weights
  # (here, all plots have a weight equal to 1)
  MC <- MetaCommunity(BCI.df)
}
```
