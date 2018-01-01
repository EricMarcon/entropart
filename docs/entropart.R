## ----knitr, echo=FALSE, cache=FALSE-----------------------
library("knitr")
opts_knit$set(concordance=TRUE)
# set global chunk options
opts_chunk$set(cache=TRUE, warning=FALSE, tidy=TRUE, fig.width=8, fig.height=6, out.width='.6\\maxwidth', tidy.opts=list(blank=FALSE, width.cutoff=60), size="scriptsize")
options(width=60)
par(mar=c(0,0,0,0))

## ----MetaCommunity----------------------------------------
library("entropart")
df <- data.frame(C1=c(10, 10, 10, 10), C2=c(0, 20, 35, 5), C3=c(25, 15, 0, 2))
row.names(df) <- c("sp1", "sp2", "sp3", "sp4")
df
w <- c(1, 2, 1)
MC <- MetaCommunity(Abundances=df, Weights=w)

## ----MCPlot, echo=FALSE, results='hide'-------------------
par(mar=c(2.5, 4, 1, 0))
plot(MC)

## ----Ps---------------------------------------------------
MC$Ps

## ----data-------------------------------------------------
data("Paracou618")
summary(Paracou618.MC)

## ----AbdVector--------------------------------------------
data("Paracou618")
PAbd <- as.AbdVector(Paracou618.MC$Ns)

## ----ProbaVector------------------------------------------
PProba <- as.ProbaVector(Paracou618.MC$Ps)

## ----ProbaVectorPlot, echo=FALSE, results='hide'----------
plot(PAbd)

## ----lnq, echo=FALSE, results='hide'----------------------
curve(log(x), 0, 1, lty=1, ylab = expression(ln[q](x)))
curve(lnq(x, 0), 0, 1, lty = 2, add = TRUE)
curve(lnq(x, 2), 0, 1, lty = 3, add = TRUE)
curve(lnq(x, 3), 0, 1, lty = 4, add = TRUE)  
legend("bottomright", legend = c(expression(ln[0](x)), "ln(x)", expression(ln[2](x)), expression(ln[3](x))), lty = c(2, 1, 3, 4), inset=  0.02)

## ----ShannonP6.1------------------------------------------
Tsallis(PProba, q = 1)

## ----ShannonP6.2------------------------------------------
Shannon(PProba)

## ----Coverage---------------------------------------------
Coverage(PAbd)

## ----bcShannonP6.1----------------------------------------
bcTsallis(PAbd, q = 1)

## ----bcShannonP6.3----------------------------------------
Tsallis(PAbd, q = 1)

## ----bcShannonP6.4----------------------------------------
Tsallis(PProba, q = 1)

## ----Diversity--------------------------------------------
expq(Tsallis(PAbd, q = 2), q = 2)
Diversity(PAbd, q = 2)

## ----Hurlbert---------------------------------------------
Hurlbert(PProba, k = 2)
Hurlbert(PAbd, k = 2)
HurlbertD(PAbd, k = 2)

## ----AlphaEntropy-----------------------------------------
e <- AlphaEntropy(Paracou618.MC, q = 1)
summary(e)

## ----DivPart----------------------------------------------
p <- DivPart(q = 1, MC = Paracou618.MC, Biased = FALSE)
summary(p)
p$CommunityAlphaEntropies

## ----DivPartPlot, echo=FALSE, results='hide'--------------
par(mar=c(4, 4, 2, 1))
plot(p)

## ----DivEst-----------------------------------------------
de <- DivEst(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", Simulations = 100)
summary(de)

## ----DivEstPlot, echo=FALSE, results='hide', out.width='.8\\maxwidth'----
par(mar=c(4, 4, 2, 1))
plot(de)

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

