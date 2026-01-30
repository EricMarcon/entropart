# Changelog

## entropart 1.6-16

CRAN release: 2025-02-07

### External changes

- Replaced `Geom*$default_aes` by their values for compatibility with
  ggplot2 3.6.0 (PR
  [\#3](https://github.com/EricMarcon/entropart/issues/3) by
  [@teunbrand](https://github.com/teunbrand))

## entropart 1.6-15

CRAN release: 2024-08-26

### Bug correction

- [`PhyloApply()`](https://ericmarcon.github.io/entropart/reference/PhyloApply.md)
  sometimes raised an “Error in m\[, 1\] : incorrect number of
  dimensions” due to the automatic conversion of a single-row matrix to
  a vector by [`lapply()`](https://rdrr.io/r/base/lapply.html). This
  happened when the tree’s upper slice contained a single species.

### User-visible changes

- pkgdown site is now built with the bootstrap 5 template.

## entropart 1.6-12

CRAN release: 2023-05-31

### Bug correction

- Corrected
  [`Coverage()`](https://ericmarcon.github.io/entropart/reference/Coverage.md)
  documentation.
- “Marcon” correction of Shannon’s entropy never returned Grassberger’s
  estimate. This was very rarely an issue.
- Phylogenetic trees of class `phylo` with multichotomies raised an
  error when they were preprocessed by
  [`Preprocess.Tree()`](https://ericmarcon.github.io/entropart/reference/Preprocess.Tree.md).

## entropart 1.6-11

CRAN release: 2022-07-11

### Bug correction

- Corrected NEWS file structure.

## entropart 1.6-10

CRAN release: 2022-04-02

### Improvements

- Deprecated `aes_(x=~var)` aesthetics in ggplots replaced by
  `aes(x=.data$var)`.
- entropart no longer suggests SPECIES.

## entropart 1.6-9

CRAN release: 2022-03-05

### Bug correction

- `Coverage2Size` raised an error when the distribution had no
  singletons.

### Improvements

- entropart no longer depends on ggplot2.

## entropart 1.6-8

CRAN release: 2021-06-04

### Bug correction

- An error occured when a function was called with `entopart::` prefix
  and `CheckArguments = TRUE` and entropart was loaded.

## entropart 1.6-7

CRAN release: 2021-04-16

### Bug correction

- argument `main` in
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  was ignored.

### Improvements

- Continuous integration by GitHub Actions.
- [`CommunityProfile()`](https://ericmarcon.github.io/entropart/reference/CommunityProfile.md)
  does not recenter simulated diversity values if simulated community
  size is not that of the actual community. `$mid` can store mean
  simulated values.
- CodeFactor recommendations applied.
- `col`, `pch`, `cex` and `lty` arguments in `autoplot`s.
- The first column of a dataframe passed to
  [`MetaCommunity()`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md)
  with species names may be characters instead of factors.
- tibbles are accepted by
  [`MetaCommunity()`](https://ericmarcon.github.io/entropart/reference/MetaCommunity.md).
- Empty communities are detected by `Metacommunity()` and raise an error
  with an explicit message rather than returning obscure error messages
  when used by
  [`DivPart()`](https://ericmarcon.github.io/entropart/reference/DivPart.md)
  and others.
- Suggests rmarkdown (<https://github.com/yihui/knitr/issues/1864>)

## entropart 1.6-6

CRAN release: 2020-11-03

### Bug correction

- argument `main` in
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  was ignored.

## entropart 1.6-4

CRAN release: 2020-01-22

### Improvements

- No longer imports geiger package (required by CRAN).

## entropart 1.6-3

CRAN release: 2020-01-08

### New features

- Argument `JackMax` to limit the order of the jacknife estimator in
  [`Richness()`](https://ericmarcon.github.io/entropart/reference/Richness.md),
  whatever the data.

### Bug correction

- plot.SpeciesDistribution did not fit the logseries and the broken
  stick distributions correctly.
- autoplot.DivProfile did not produce appropriate labels in the local
  communities profile.

## entropart 1.6-1

CRAN release: 2019-04-11

### New features

- Estimation of diversity at a chosen level (sample size or coverage).
  `DivAccum()` function.
- Entropy accumulation functions.
- ggplot2 supported.
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods added for entropart objects.
- The “Best” estimator of diversity is now “UnveilJ” and the default
  estimator of richness is “Jackknife”.
- The “ChaoWangJost” estimator is renamed “ChaoJost”.

### Improvements

- Unit tests added.
- Vignette by pkgdown.

### Bug correction

- The jaccknife estimator of richness returned an error for communities
  where all species had the same abundance.
- `Richness` returned 0 instead of 1 for a community with a single
  species.

## entropart 1.5-3

CRAN release: 2018-02-05

### Improvements

- On Travis now.
- Reduced package size.
- The rule to calculate the number of individuals of MetaCommunities has
  been changed to improve gamma diversity bias correction. See the user
  manual vignette.
- Generic function arguments cleaned up.

### Bug Correction

- Very large metacommunities returned an integer overflow error.

## entropart 1.4-8

CRAN release: 2017-11-22

### Bug Correction

- [`HqzBeta()`](https://ericmarcon.github.io/entropart/reference/HqzBeta.md)
  returned erroneous values if a species probability was equal to zero.

### Improvements

- On GitHub now.
- Documentation updated: phylogenetic dendrograms can be of class
  `phylo`, `phylog`, `hclust` or `PPtree` whatever the function.
- The introduction vignette is HTML now.
- A new vignette is dedicated to phylogenies.

## entropart 1.4-7

CRAN release: 2017-03-29

### Bug Correction

- Argument checking (`CheckArguments = TRUE`) is not possible when the
  package is not loaded and a function is called by
  `entropart::function()`. An error was returned. It is replaced by a
  warning.

### Improvements

- Explicit export of all non-internal functions instead of
  `exportPattern("^[[:alpha:]]+")`
- Updated references to published articles.
- Updated
  [`help("entropart")`](https://ericmarcon.github.io/entropart/reference/entropart-package.md).
- New introduction vignette.
- Vignettes compiled with *knitr* instead of *Sweave*.

## entropart 1.4-6

CRAN release: 2016-12-09

### Improvements

- LazyData is used to save memory.
- Better reporting of the argument names in embedded calls of functions.

### Bug Correction

- The simulation of log-series communities was incorrect.

## entropart 1.4-5

CRAN release: 2016-02-17

### User-visible changes

- Generalized Simpson’s entropy and diversity added (`GenSimpson` and
  `GenSimpsonD`).
- `Originality.Species()` is deprecated because it is pointless.
  [`ade4::originality()`](https://adeverse.github.io/ade4/reference/originality.html)
  allows calculating it for q=2. Leinster (2009) and Leinster and
  Meckes (2015) showed that `Originality.Species()` does not depend on
  the order of diversity.

### Improvements

- ZhangGrabchak estimator of entropy is now calculated by the C code of
  [`EntropyEstimation::Tsallis.z`](https://rdrr.io/pkg/EntropyEstimation/man/Tsallis.z.html)/`Entropy.z`
  rather than the R code of
  [`bcTsallis()`](https://ericmarcon.github.io/entropart/reference/Tsallis.md).
  This is much faster when the number of individual is high. Applies to
  `ChaoWangJost` (Best) estimator too.

## entropart 1.4-4

CRAN release: 2016-01-13

### User-visible changes

- [`DivProfile()`](https://ericmarcon.github.io/entropart/reference/DivProfile.md)
  now allows computing bootstrap confidence intervals.

### Bug Correction

- The entropy estimation (of order different from 1) of a distribution
  with no singleton returned `NA` with `ChaoWangJost` correction.
  Reported by Zach Marion. Only partly corrected in Version 1.4-1.
  Corrected.
- `DivEst` returned incorrect beta diversity if q was not 1. Corrected.

## entropart 1.4-3

CRAN release: 2015-12-07

### User-visible changes

- All scalar values of diversity or entropy are now named. Their name is
  the bias correction used to obtain them.
- The `Unveiled` estimator is more versatile. `Correction = "Unveil"` is
  deprecated and replaced by `UnveilC`, `UnveiliC` or `UnveilJ` in
  functions such as
  [`Tsallis()`](https://ericmarcon.github.io/entropart/reference/Tsallis.md)
  or
  [`Diversity()`](https://ericmarcon.github.io/entropart/reference/Diversity.md).

### Improvements

- Parallelization of
  [`DivProfile()`](https://ericmarcon.github.io/entropart/reference/DivProfile.md),
  [`CommunityProfile()`](https://ericmarcon.github.io/entropart/reference/CommunityProfile.md)
  and
  [`PhyloApply()`](https://ericmarcon.github.io/entropart/reference/PhyloApply.md)
  using the parallel package *mclapply*. No effect on Windows, pretty
  much faster on other systems.
- Extensive use of [`vapply()`](https://rdrr.io/r/base/lapply.html)
  instead of [`sapply()`](https://rdrr.io/r/base/lapply.html) makes some
  functions faster.
- [`AllenH()`](https://ericmarcon.github.io/entropart/reference/AllenH.md)
  and
  [`ChaoPD()`](https://ericmarcon.github.io/entropart/reference/ChaoPD.md)
  returned `NA` if the tree contained more species than the probability
  vector. Now, the tree may be pruned or kept unchanged and extra
  species considered to have probabilities 0.

### Bug Correction

- Using `phylog` trees in `AllenH` and
  [`ChaoPD()`](https://ericmarcon.github.io/entropart/reference/ChaoPD.md)
  returned erroneous unnormalized diversity (divided by two) because of
  the conVersion of `phylog` to `htree` divides branch lengths by two.
  Corrected.
- The richness estimator `iChao1` returned `NA` if the distibution
  contained singletons but no doubletons. Corrected.

## entropart 1.4.1

CRAN release: 2015-09-29

### New Features

- `phylog` objects (deprecated in *ade4*) are replaced by `phylo` trees
  from package *ape* in the definition of the `PPtree` class. Issues
  caused by `phylog` such as replacing `.` and `-` by `_` in species
  names do not occur any longer. `phylog` trees are still accepted for
  compatibility.
- [`ChaoPD()`](https://ericmarcon.github.io/entropart/reference/ChaoPD.md)
  and
  [`AllenH()`](https://ericmarcon.github.io/entropart/reference/AllenH.md)
  now accept `phylo` trees.
- `Richness` now returns a named value. The name contains the estimator
  used.
- Updated *CITATION*: the paper about this package has been published:
  Eric Marcon, Bruno Herault (2015). entropart: An R Package to Measure
  and Partition Diversity. *Journal of Statistical Software*, 67(8),
  1-26.

### Bug Correction

- The entropy estimation of a distribution with no singleton returned
  `NA` with `ChaoWangJost` correction. Corrected.
- Entropy or diversity of a vector of zeros returned 0. It now returns
  `NA`.

## entropart 1.3.3

CRAN release: 2015-09-08

### New Features

- Abundance and probability vector objects. See
  [`?SpeciesDistribution`](https://ericmarcon.github.io/entropart/reference/SpeciesDistribution.md).
- Hurlbert diversity. See
  [`?Hurlbert`](https://ericmarcon.github.io/entropart/reference/Hurlbert.md).
- `Optimal.Similarity`.
- Miller-Madow estimator of entropy (Miller, 1955) added in
  [`bcShannon()`](https://ericmarcon.github.io/entropart/reference/Shannon.md).
- Chao and Jost (2015) estimator of diversity added in
  [`bcTsallis()`](https://ericmarcon.github.io/entropart/reference/Tsallis.md)
  and
  [`bcDiversity()`](https://ericmarcon.github.io/entropart/reference/Diversity.md).
  New “best” estimator.
- Chao et al. (2015) probability estimation of observed species. See
  `?TunedPs`.
- Estimators of the number of species. See
  [`?Richness`](https://ericmarcon.github.io/entropart/reference/Richness.md).
- Abundance Frequency Count of species. See
  [`?AbdFreqCount`](https://ericmarcon.github.io/entropart/reference/AbdFreqCount.md).
- Community profiles can be calculated with confidence intervals. See
  [`?CommunityProfile`](https://ericmarcon.github.io/entropart/reference/CommunityProfile.md).
- Random Communities. See
  [`?rCommunity`](https://ericmarcon.github.io/entropart/reference/rCommunity.md).

### Bug Correction

- Applying `bcTsallis` and similar functions with a probability vector
  instead of abundance values could cause errors depending in the
  correction. Correction is now forced to `None` with a warning.
- Allowed rounding error was too small on some systems (typically
  r-patched-solaris-sparc) to recognize probability vectors. The
  difference between their sum and 1 had to be less than 3 times
  `.Machine$double.eps`. Now set to S times (where S is the number of
  species, i.e the vector’s length).

## entropart 1.2.1

CRAN release: 2014-11-13

### New Features

- Zhang and Grabchak (2014) bias correction for Shannon beta entropy
  added.
- Unbiased estimator of Rao’s entropy added (`bcRao`).

### Bug Correction

- `DqZ()` and
  [`Hqz()`](https://ericmarcon.github.io/entropart/reference/Hqz.md)
  returned an error if all probability values were 0 except one.

### Improvements

- Improved readability of error messages for bad arguments.
- Improved formating of `summmary.DivPart()`. Lines were too long.
- Improved legend for the x-axis of `plot.DivPart` (“alpha and gamma”
  instead of “alpha/gamma”).
- Improved support of `PhyloValue` objects (summary added).
- Improved help for `MetaCommunity`.

## entropart 1.2.0

CRAN release: 2014-05-17

### Bug Correction

- [`ChaoPD()`](https://ericmarcon.github.io/entropart/reference/ChaoPD.md)
  returned an incorrect value when q=0 and some probabilities =0.

## entropart 1.2.0

CRAN release: 2014-05-17

### New Features

- Full support of similarity-based diversity added

### Improvements

- Default values for arguments added whenever possible.

## entropart 1.1.4

### New Features

- Zhang(2012) bias correction for Shannon entropy added.
- Zhang and Grabchak (2014) bias correction for Tsallis entropy added.

### Bug Correction

- `Divest()` always calculated neutral diversity of simulated
  communities so the confidence interval was erroneous for
  phylodiversity. Corrected.

## entropart 1.1.3

CRAN release: 2014-03-22

### New Features

- `Paracou618.dist` distance matrix between species of `Paracou618.MC`
  added.
- Chao, Wang and Jost (2013) bias correction for Shannon entropy added.
- `EntropyCI` function added: Entropy of Monte-Carlo simulated
  communities.
- Tools to manipulate MetaCommunity objects added (see
  [`?MergeMC`](https://ericmarcon.github.io/entropart/reference/MergeMC.md)).
- `SimTest` class added to test a value against a simulated distribution
  (see
  [`?SimTest`](https://ericmarcon.github.io/entropart/reference/SimTest.md)).
- Vignette added.

### Improvements

- `Imports` directive rather than `Depends` for *ade4*.
- `mergeandlabel` does not return warnings any longer (column names are
  better addressed).

### Bug Correction

- [`Hqz()`](https://ericmarcon.github.io/entropart/reference/Hqz.md) was
  erroneous for q\<\>1. Corrected.
- [`bcPhyloEntropy()`](https://ericmarcon.github.io/entropart/reference/PhyloEntropy.md)
  and
  [`bcPhyloDiversity()`](https://ericmarcon.github.io/entropart/reference/PhyloDiversity.md)
  returned an incorrect `$Distribution` component. Corrected.
- [`summary.MCentropy()`](https://ericmarcon.github.io/entropart/reference/MCentropy.md)
  did not return the name of the tree. Corrected.
- Legend was not displayed in
  `plot.DivProfile(..., Which="Communities")`. Corrected.

## entropart 1.1.1

CRAN release: 2013-09-10

- First Version.
