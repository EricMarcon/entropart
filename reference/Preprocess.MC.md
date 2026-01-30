# Preprocessing of metacommunity data

Calculates summary statistics of a metacommunity

## Usage

``` r
Preprocess.MC(Nsi, Wi)
```

## Arguments

- Nsi:

  A matrix containing abundance data, species in line, communities in
  column.

- Wi:

  Normalized community weights (sum to 1).

## Details

This is an internal function. It is used by
[`MetaCommunity`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md)
for real data and
[`DivEst`](https://ericmarcon.github.io/entropart/reference/DivEst.md)
for simulations.

## Value

An object of class **MetaCommunity**. It is a list:

- Nsi:

  A matrix containing abundance data, species in line, communities in
  column

- Nspecies:

  The number of species

- Ncommunities:

  The number of communities

- Wi:

  Normalized community weights (sum to 1)

- N:

  The total number of individuals

- Ni:

  A vector containing the number of individuals of each community

- Ns:

  A vector containing the number of individuals of each species

- Psi:

  A matrix whose columns are the probability vectors of communities
  (each of them sums to 1)

- Ps:

  A vector containing the probability vector of the metacommunity

- SampleCoverage:

  The sample coverage of the metacommunity

- SampleCoverage.communities:

  A vector containing the sample coverages of each community

## Author

Eric Marcon \<Eric.Marcon@ecofog.gf\>
