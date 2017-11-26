# Entropy Partitioning to Measure Diversity

[![CRAN version](http://www.r-pkg.org/badges/version/entropart)](https://cran.r-project.org/package=entropart)
[![](http://cranlogs.r-pkg.org/badges/entropart)](https://cran.r-project.org/package=entropart)
[![Research software impact](http://depsy.org/api/package/cran/entropart/badge.svg)](http://depsy.org/package/r/entropart)

entropart is an R package that provides functions to calculate alpha, beta and gamma diversity of communities, 
including phylogenetic and functional diversity.
  
Estimation-bias corrections are available.

# Details

In the entropart package, individuals of different *species* are counted in several *communities* which may (or not) 
be agregated to define a *metacommunity*. 
In the metacommunity, the probability to find a species in the weighted average of probabilities in communities. 
This is a naming convention, which may correspond to plots in a forest inventory or any data organized the same way.

Basic functions allow computing diversity of a community. 
Data is simply a vector of probabilities (summing up to 1) or of abundances (integer values that are numbers of individuals). 
Calculate entropy with functions such as *Tsallis*, *Shannon*, *Simpson*, *Hurlbert* or *GenSimpson* 
and explicit diversity (i.e. effective number of species) with *Diversity* and others. 
By default, the best available estimator of diversity will be used, according to the data.
  
Communities can be simulated by *rCommunity*, explicitely declared as a species distribution (*as.AbdVector* or *as.ProbaVector*), 
and plotted.
  
Phylogenetic entropy and diversity can be calculated if a phylogenetic (or functional), ultrametric tree is provided. 
See *PhyloEntropy*, *Rao* for examples of entropy and *PhyloDiversity* to calculate phylodiversity, 
with the state-of-the-art estimation-bias correction. 
Similarity-based diversity is calculated with *Dqz*, based on a similarity matrix.

A full documentation is available in the vignette. 
Type: *vignette("entropart")*. 
A quick introuction is in *vignette("Introduction", "entropart")*.

## Reference

Marcon, E. and Herault, B. (2015). entropart: An R Package to Measure and Partition Diversity.
*Journal of Statistical Software*. 67(8): 1-26.
