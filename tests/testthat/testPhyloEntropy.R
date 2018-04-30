testthat::context("PhyloEntropy")

# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Species probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ns)
# Taxonomy: Build a ppTree
ppTree <- Preprocess.Tree(Paracou618.Taxonomy)
# Build a phylog
phylogTree <- ade4::hclust2phylog(ppTree$hTree)


# Check PhyloEntropy does not change with dendrogram format
testthat::test_that("PhyloEntropy does not depend on tree format", {
  testthat::skip_on_cran()
  # phylog and phylo
  testthat::expect_equal(AlphaEntropy(Paracou618.MC, q=1, Tree = phylogTree)$Total,
               AlphaEntropy(Paracou618.MC, q=1, Tree = Paracou618.Taxonomy)$Total)
})


# Check PhyloEntropy of order 2 equals Rao
testthat::test_that("PhyloEntropy of order 2 equals Rao", {
  testthat::skip_on_cran()
  # No correction
  testthat::expect_equal(as.numeric(PhyloEntropy(Ps, 2, Paracou618.Taxonomy, Normalize = FALSE)$Total), 
               as.numeric(Rao(Ps, Paracou618.Taxonomy)))
  # With correction
  testthat::expect_equal(as.numeric(PhyloEntropy(Ns, 2, Paracou618.Taxonomy, Normalize = FALSE, Correction = "ChaoJost")$Total), 
               as.numeric(Rao(Ns, Paracou618.Taxonomy)))
  # HqZ
  testthat::expect_equal(as.numeric(PhyloEntropy(Ps, 2, Paracou618.Taxonomy, Normalize = TRUE)$Total), 
               as.numeric(Hqz(Ps, 2, Z=1-as.matrix(phylogTree$Wdist^2/6))))
})


# Check PhyloDiversity equals ChaoPD
testthat::test_that("PhyloDiversity equals ChaoPD", {
  testthat::skip_on_cran()
  # Order 2
  testthat::expect_equal(as.numeric(PhyloDiversity(Ps, 2, Paracou618.Taxonomy)$Total), 
               as.numeric(ChaoPD(Ps, 2, Paracou618.Taxonomy)))
  # Order 1
  testthat::expect_equal(as.numeric(PhyloDiversity(Ps, 1, Paracou618.Taxonomy)$Total), 
               as.numeric(ChaoPD(Ps, 1, Paracou618.Taxonomy)))
})


# Check PhyloEntropy of order 2 of a star dendrogram equals Simpson
# Create a star dendrogram (actually, almost a star: phylo trees must be binary)
NewickStar <- "(A:1,(B:.99999,C:.99999):0.00001):0;"
phyStar <- ape::read.tree(text=NewickStar)
# Randomly assign probabilities to 3 species
BrockenStick <- sort(runif(2))
PsStar <- c(BrockenStick[1], BrockenStick[2]-BrockenStick[1], 1-BrockenStick[2])
names(PsStar) <- c("A", "B", "C")

testthat::test_that("PhyloEntropy of order 2 of a star dendrogram equals Simpson", {
  testthat::skip_on_cran()
  # PhyloEntropy vs Rao
  testthat::expect_equal(as.numeric(PhyloEntropy(PsStar, q=2, Tree=phyStar, Normalize = FALSE)$Total), 
               as.numeric(Rao(PsStar, phyStar)), tolerance = 100*sqrt(.Machine$double.eps))
  # Rao vs Simpson. Should be equal but the tree is not exactly a star.
  testthat::expect_lt((Simpson(PsStar)-Rao(PsStar, phyStar))/Simpson(PsStar), 
             1/1000)
})

