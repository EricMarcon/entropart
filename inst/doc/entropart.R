## ----knitr, echo=FALSE, cache=FALSE-----------------------
library("knitr")
opts_knit$set(concordance=TRUE)
# set global chunk options
opts_chunk$set(cache=TRUE, warning=FALSE, tidy=TRUE, fig.width=8, fig.height=6, out.width='.4\\maxwidth', tidy.opts=list(blank=FALSE, width.cutoff=60), size="scriptsize")
options(width=60)
par(mar=c(0,0,0,0))

## ----DivProfile-------------------------------------------
dp <- DivProfile(seq(0, 2, 0.2), Paracou618.MC, Biased = FALSE, NumberOfSimulations = 10)
summary(dp)

## ----DivProfilePlot, echo=FALSE, results='hide', out.width='.8\\maxwidth'----
plot(dp)

## ----ShannonBeta------------------------------------------
ShannonBeta(Paracou618.MC$Psi[, 1], PProba)

## ----PhyloDiversity---------------------------------------
phd <- bcPhyloDiversity(PAbd, q = 1, Tree = Paracou618.Taxonomy, Normalize = TRUE)
summary(phd)

## ----PhyloDiversityPlot, echo=FALSE, results='hide'-------
par(mar=c(4, 4, 2, 1))
plot(phd, main = "")

## ----PhylodivPart-----------------------------------------
dp <- DivPart(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", Tree = Paracou618.Taxonomy)
summary(dp)

## ----PhyloBetaEntropy-------------------------------------
summary(BetaEntropy(Paracou618.MC, q = 2, Tree = Paracou618.Taxonomy, Correction = "None", Normalize = FALSE))

## ----divc-------------------------------------------------
library("ade4")
divc(as.data.frame(Paracou618.MC$Wi), disc(as.data.frame(Paracou618.MC$Nsi), Paracou618.Taxonomy$Wdist))

## ----Dqz--------------------------------------------------
DistanceMatrix <- as.matrix(Paracou618.dist)
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
bcDqz(PAbd, q = 2, Z)

## ----DqzHCDT----------------------------------------------
Dqz(PProba, q = 2, Z = diag(length(PProba)))
Diversity(PProba, q = 2)

## ----Hqz--------------------------------------------------
Hqz(PProba, q = 2, Z)
lnq(Dqz(PProba, q = 2, Z), q = 2)

## ----AlphaEntropyZ----------------------------------------
e <- AlphaEntropy(Paracou618.MC, q = 1, Z = Z)
summary(e)

## ----rCommunity, fig.show='hide'--------------------------
rCommunity(n = 1, size = 1000, S=300, Distribution = "lnorm", sd=1) -> NsRef
Richness(as.ProbaVector(NsRef))
Richness(NsRef)
plot(NsRef, Distribution="lnorm")

## ----rCommunityFig, echo=FALSE, results='hide'------------
plot(NsRef, Distribution="lnorm") -> fit

## ----EntropyCI--------------------------------------------
SimulatedDiversity <- expq(EntropyCI(FUN = Tsallis, Simulations = 100, Ns = PAbd, q = 1), q = 1)
Diversity(PAbd, q = 1)
quantile(SimulatedDiversity, probs = c(0.025, 0.975))

## ----CommunityProfile-------------------------------------
bcProfile <- CommunityProfile(Diversity, PAbd, NumberOfSimulations = 10)
Profile <- CommunityProfile(Diversity, PProba)

## ----CommunityProfileFigCode, eval=FALSE------------------
#  plot(bcProfile)
#  CEnvelope(Profile, lty=3)
#  legend("topright", c("Bias Corrected", "Biased"), lty=c(1,3), inset=0.02)

## ----CommunityProfileFig, echo=FALSE, results='hide'------
par(mar=c(4, 4, 2, 1))
plot(bcProfile)
CEnvelope(Profile, lty=3)
legend("topright", c("Bias Corrected", "Biased"), lty=c(1,3), inset=0.02)

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

