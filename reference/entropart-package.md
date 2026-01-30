# Entropy Partitioning to Measure Diversity

Functions to calculate alpha, beta and gamma diversity of communities,
including phylogenetic and functional diversity.

Estimation-bias corrections are available.

## Details

In the entropart package, individuals of different "species" are counted
in several "communities" which may (or not) be agregated to define a
"metacommunity". In the metacommunity, the probability to find a species
in the weighted average of probabilities in communities. This is a
naming convention, which may correspond to plots in a forest inventory
or any data organized the same way.

Basic functions allow computing diversity of a community. Data is simply
a vector of probabilities (summing up to 1) or of abundances (integer
values that are numbers of individuals). Calculate entropy with
functions such as
[`Tsallis`](https://ericmarcon.github.io/entropart/reference/Tsallis.md),
[`Shannon`](https://ericmarcon.github.io/entropart/reference/Shannon.md),
[`Simpson`](https://ericmarcon.github.io/entropart/reference/Simpson.md),
[`Hurlbert`](https://ericmarcon.github.io/entropart/reference/Hurlbert.md)
or
[`GenSimpson`](https://ericmarcon.github.io/entropart/reference/GenSimpson.md)
and explicit diversity (i.e. effective number of species) with
[`Diversity`](https://ericmarcon.github.io/entropart/reference/Diversity.md)
and others. By default, the best available estimator of diversity will
be used, according to the data.

Communities can be simulated by
[`rCommunity`](https://ericmarcon.github.io/entropart/reference/rCommunity.md),
explicitely declared as a species distribution
([`as.AbdVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md)
or
[`as.ProbaVector`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md)),
and plotted.

Phylogenetic entropy and diversity can be calculated if a phylogenetic
(or functional), ultrametric tree is provided. See
[`PhyloEntropy`](https://ericmarcon.github.io/entropart/reference/PhyloEntropy.md),
[`Rao`](https://ericmarcon.github.io/entropart/reference/Rao.md) for
examples of entropy and
[`PhyloDiversity`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)
to calculate phylodiversity, with the state-of-the-art estimation-bias
correction. Similarity-based diversity is calculated with
[`Dqz`](https://ericmarcon.github.io/entropart/reference/Dqz.md), based
on a similarity matrix.

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

The last line of the code calls the
[`MetaCommunity`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md)
function to create an object that will be used by all metacommunity
functions, such as
[`DivPart`](https://ericmarcon.github.io/entropart/reference/DivPart.md)
(to partition diversity),
[`DivEst`](https://ericmarcon.github.io/entropart/reference/DivEst.md)
(to partition diversity and calculate confidence interval of its
estimation) or
[`DivProfile`](https://ericmarcon.github.io/entropart/reference/DivProfile.md)
(to compute diversity profiles).

A full documentation is available in the vignette. Type:
[`vignette("entropart")`](https://ericmarcon.github.io/entropart/articles/entropart.md).
A quick introuction is in `vignette("introduction", "entropart")`.

## Author

Eric Marcon, Bruno Herault

## References

Grabchak, M., Marcon, E., Lang, G., and Zhang, Z. (2017). The
Generalized Simpson's Entropy is a Measure of Biodiversity. *Plos One*,
12(3): e0173305.

Marcon, E. (2015) Practical Estimation of Diversity from Abundance Data.
*HAL* 01212435: 1-27.

Marcon, E. and Herault, B. (2015). entropart: An R Package to Measure
and Partition Diversity. *Journal of Statistical Software*, 67(8): 1-26.

Marcon, E., Herault, B. (2015). Decomposing Phylodiversity. *Methods in
Ecology and Evolution* 6(3): 333-339.

Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The
Decomposition of Shannon's Entropy and a Confidence Interval for Beta
Diversity. *Oikos* 121(4): 516-522.

Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang G. (2014).
Generalization of the partitioning of Shannon diversity. *PLOS One*
9(3): e90289.

Marcon, E., Zhang, Z. and Herault, B. (2014). The decomposition of
similarity-based diversity and its bias correction. *HAL*
hal-00989454(version 3).
