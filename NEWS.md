# entropart 1.5-3-9008

## New features

- Estimation of diversity at a chosen level (sample size or coverage).
- Entropy accumulation functions.
- ggplot2 supported. `autoplot()` methods added for entropart objects.
- The "Best" estimator of diversity is now "UnveilJ" and the default estimator of richness is "Jackknife".
- The "ChaoWangJost" estimator is renamed "ChaoJost".

## Improvements

- Unit tests added.
- Vignette by pkgdown.



# entropart 1.5-3

## Improvements

- On Travis now.
- Reduced package size.
- The rule to calculate the number of individuals of MetaCommunities has been changed to improve gamma diversity bias correction. See the user manual vignette.
- Generic function arguments cleaned up.

## Bug Correction

- Very large metacommunities returned an integer overflow error.



# entropart 1.4-8

## Bug Correction

- `HqzBeta()` returned erroneous values if a species probability was equal to zero.

## Improvements

- On GitHub now.
- Documentation updated: phylogenetic dendrograms can be of class `phylo`, `phylog`, `hclust` or `PPtree` whatever the function.
- The introduction vignette is HTML now.
- A new vignette is dedicated to phylogenies.



# entropart 1.4-7

## Bug Correction

- Argument checking (`CheckArguments = TRUE`) is not possible when the package is not loaded and a function is called by `entropart::function()`. An error was returned. It is replaced by a warning.
  
## Improvements

- Explicit export of all non-internal functions instead of `exportPattern("^[[:alpha:]]+")`
- Updated references to published articles.
- Updated `help("entropart")`.
- New introduction vignette.
- Vignettes compiled with _knitr_ instead of _Sweave_.
  


# entropart 1.4-6

# Improvements

- LazyData is used to save memory.
- Better reporting of the argument names in embedded calls of functions.

## Bug Correction

- The simulation of log-series communities was incorrect.



# entropart 1.4-5

## User-visible changes

- Generalized Simpson's entropy and diversity added (`GenSimpson` and `GenSimpsonD`).
- `Originality.Species()` is deprecated because it is pointless. `ade4::originality()` allows calculating it for q=2. Leinster (2009) and Leinster and Meckes (2015) showed that `Originality.Species()` does not depend on the order of diversity.
  
## Improvements

- ZhangGrabchak estimator of entropy is now calculated by the C code of `EntropyEstimation::Tsallis.z`/`Entropy.z` rather than the R code of `bcTsallis()`. This is much faster when the number of individual is high. Applies to `ChaoWangJost` (Best) estimator too.



# entropart 1.4-4

## User-visible changes 

- `DivProfile()` now allows computing bootstrap confidence intervals.


## Bug Correction

- The entropy estimation (of order different from 1) of a distribution with no singleton returned `NA` with `ChaoWangJost` correction. Reported by Zach Marion. Only partly corrected in Version 1.4-1. Corrected.
- `DivEst` returned incorrect beta diversity if q was not 1. Corrected.



# entropart 1.4-3

## User-visible changes 

- All scalar values of diversity or entropy are now named. Their name is the bias correction used to obtain them.
- The `Unveiled` estimator is more versatile. `Correction = "Unveil"` is deprecated and replaced by `UnveilC`, `UnveiliC` or `UnveilJ` in functions such as `Tsallis()` or `Diversity()`.

# Improvements 

- Parallelization of `DivProfile()`, `CommunityProfile()` and `PhyloApply()` using the parallel package _mclapply_. No effect on Windows, pretty much faster on other systems.
- Extensive use of `vapply()` instead of `sapply()` makes some functions faster.
- `AllenH()` and `ChaoPD()` returned `NA` if the tree contained more species than the probability vector. Now, the tree may be pruned or kept unchanged and extra species considered to have probabilities 0.

## Bug Correction

- Using `phylog` trees in `AllenH` and `ChaoPD()` returned erroneous unnormalized diversity (divided by two) because of the conVersion of `phylog` to `htree` divides branch lengths by two. Corrected.
- The richness estimator `iChao1` returned `NA` if the distibution contained singletons but no doubletons. Corrected.




# entropart 1.4.1


## New Features

- `phylog` objects (deprecated in _ade4_) are replaced by `phylo` trees from package _ape_ in the definition of the `PPtree` class. Issues caused by `phylog` such as replacing `.` and `-` by `_` in species names do not occur any longer. `phylog` trees are still accepted for compatibility.
- `ChaoPD()` and `AllenH()` now accept `phylo` trees.
- `Richness` now returns a named value. The name contains the estimator used.
- Updated _CITATION_: the paper about this package has been published: Eric Marcon, Bruno Herault (2015). entropart: An R Package to Measure and Partition Diversity. _Journal of Statistical Software_, 67(8), 1-26.
  
## Bug Correction 

- The entropy estimation of a distribution with no singleton returned `NA` with `ChaoWangJost` correction. Corrected.
- Entropy or diversity of a vector of zeros returned 0. It now returns `NA`.
  
  
  
# entropart 1.3.3

## New Features

- Abundance and probability vector objects. See `?SpeciesDistribution`.
- Hurlbert diversity. See `?Hurlbert`.
- `Optimal.Similarity`.
- Miller-Madow estimator of entropy (Miller, 1955) added in `bcShannon()`.
- Chao and Jost (2015) estimator of diversity added in `bcTsallis()` and `bcDiversity()`. New "best" estimator.
- Chao et al. (2015) probability estimation of observed species. See `?TunedPs`.
- Estimators of the number of species. See `?Richness`.
- Abundance Frequency Count of species. See `?AbdFreqCount`.
- Community profiles can be calculated with confidence intervals. See `?CommunityProfile`.
- Random Communities. See `?rCommunity`.

## Bug Correction 

- Applying `bcTsallis` and similar functions with a probability vector instead of abundance values could cause errors depending in the correction. Correction is now forced to `None` with a warning.
- Allowed rounding error was too small on some systems (typically r-patched-solaris-sparc) to recognize probability vectors. The difference between their sum and 1 had to be less than 3 times `.Machine$double.eps`. Now set to S times (where S is the number of species, i.e the vector's length).



# entropart 1.2.1

## New Features

- Zhang and Grabchak (2014) bias correction for Shannon beta entropy added.
- Unbiased estimator of Rao's entropy added (`bcRao`).

## Bug Correction 

- `DqZ()` and `Hqz()` returned an error if all probability values were 0 except one.

## Improvements

- Improved readability of error messages for bad arguments.
- Improved formating of `summmary.DivPart()`. Lines were too long.
- Improved legend for the x-axis of `plot.DivPart` ("alpha and gamma" instead of "alpha/gamma").
- Improved support of `PhyloValue` objects (summary added).
- Improved help for `MetaCommunity`.



# entropart 1.2.0

## Bug Correction 

- `ChaoPD()` returned an incorrect value when q=0 and some probabilities =0.



# entropart 1.2.0

## New Features

- Full support of similarity-based diversity added 

## Improvements

- Default values for arguments added whenever possible.


# entropart 1.1.4

## New Features

- Zhang(2012) bias correction for Shannon entropy added. 
- Zhang and Grabchak (2014) bias correction for Tsallis entropy added. 

## Bug Correction 

- `Divest()` always calculated neutral diversity of simulated communities so the confidence interval was erroneous for phylodiversity. Corrected.



# entropart 1.1.3

## New Features

- `Paracou618.dist` distance matrix between species of `Paracou618.MC` added. 

# Improvements

- `Imports` directive rather than `Depends` for _ade4_.
- `mergeandlabel` does not return warnings any longer (column names are better addressed).

## Bug Correction

- Legend was not displayed in `plot.DivProfile(..., Which="Communities")`. Corrected.



# entropart 1.1.2

## New Features
  
- Chao, Wang and Jost (2013) bias correction for Shannon entropy added.
- `EntropyCI` function added: Entropy of Monte-Carlo simulated communities.
- Tools to manipulate MetaCommunity objects added (see `?MergeMC`).
- `SimTest` class added to test a value against a simulated distribution (see `?SimTest`).
- Vignette added.

## Bug Correction 

- `Hqz()` was erroneous for q<>1. Corrected.
- `bcPhyloEntropy()` and `bcPhyloDiversity()` returned an incorrect `$Distribution` component. Corrected.
- `summary.MCentropy()` did not return the name of the tree. Corrected.


# entropart 1.1.1
  
- First Version.
