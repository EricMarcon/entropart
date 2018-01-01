## ----knitr, echo=FALSE, cache=FALSE-----------------------
library("knitr")
opts_knit$set(concordance=TRUE)
# set global chunk options
opts_chunk$set(cache=TRUE, warning=FALSE, tidy=TRUE, fig.width=8, fig.height=6, out.width='.6\\maxwidth', tidy.opts=list(blank=FALSE, width.cutoff=60), size="scriptsize")
options(width=60)
par(mar=c(0,0,0,0))

## ----PhyloApply-------------------------------------------
pa <- PhyloApply(Tree=Paracou618.Taxonomy, FUN=bcTsallis, NorP=PAbd)
summary(pa)
exp(pa$Cuts)
exp(pa$Total)

## ----MC1--------------------------------------------------
(df <- data.frame(C1 = c(10, 10, 10, 10), C2 = c(0, 20, 35, 5),
  C3 = c(25, 15, 0, 2), row.names = c("sp1", "sp2", "sp3", "sp4")))
w <- colSums(df)
MC1 <- MetaCommunity(Abundances = df, Weights = w)

## ----MC2--------------------------------------------------
(df <- data.frame(C1 = c(10, 4), C2 = c(3, 4), row.names = c("sp1", "sp5")))
w <- colSums(df)
MC2 <- MetaCommunity(Abundances = df, Weights = w)

## ----MCMergeC---------------------------------------------
mergedMC1 <- MergeC(list(MC1, MC2))
mergedMC1$Nsi

## ----MCMergeMC--------------------------------------------
mergedMC2 <- MergeMC(list(MC1, MC2), Weights = sapply(list(MC1, MC2), function(x) (x$N)))
mergedMC2$Nsi

## ----MCMergeHi--------------------------------------------
dpAll <- DivPart(q=1, MC=mergedMC2)
summary(dpAll)

## ----MCMergeLo1-------------------------------------------
dpMC1 <- DivPart(q=1, MC=MC1)
summary(dpMC1)

