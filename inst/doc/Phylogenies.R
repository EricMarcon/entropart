## ----DistMatrix------------------------------------------------------------
dSp <- matrix(c(0, 1, 2, 1, 0, 2, 2, 2, 0), nrow=3, byrow=TRUE)
row.names(dSp) <- colnames(dSp) <- c("A", "B", "C")
dSp

## ----hclust----------------------------------------------------------------
require("stats")
plot(hTree <- hclust(as.dist(dSp), method="average"), hang=-0.01, axes = F)
axis(2)

## ----hclust2---------------------------------------------------------------

hTree$height

## ----phylo-----------------------------------------------------------------
require("ape")
plot(phyloTree <- as.phylo(hTree))
axis(1)

## ----phylo2----------------------------------------------------------------
phyloTree$edge.length

## ----phylo3----------------------------------------------------------------
phyloTree$edge.length <- 2*phyloTree$edge.length
plot(phyloTree)
axis(1)

## ----phylog----------------------------------------------------------------
require("ade4")
plot(phylogTree <- hclust2phylog(hTree))
axis(1)

## ----phylog2---------------------------------------------------------------
phylogTree$droot
phylogTree$Wdist^2/2

## ----Newick----------------------------------------------------------------
NewickABC <- "(C:2,(A:1,B:1):1);"
plot(phyloABC <- ape::read.tree(text=NewickABC))
axis(1)

