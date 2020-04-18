library("entropart")
read.csv2("L4C2.csv") -> dfMC
q <- 1

plot(MC <- MetaCommunity(dfMC))
plot(dp <- DivPart(q, MC))
Alpha <- dp$TotalAlphaDiversity
Gamma <- dp$GammaDiversity
Beta <- dp$TotalBetaDiversity

# Espace
plot(MCS1 <- MetaCommunity(dfMC[, 1:2]))
plot(dpS1 <- DivPart(q, MCS1))
BetaS1 <- dpS1$GammaDiversity/dpS1$TotalAlphaDiversity
plot(MCS2 <- MetaCommunity(dfMC[, 3:4]))
plot(dpS2 <- DivPart(q, MCS2))
BetaS2 <- dpS2$GammaDiversity/dpS2$TotalAlphaDiversity
BetaST <- sqrt(BetaS1*BetaS2)
plot(MergeMC(list(MCS1, MCS2)) -> MCS)
plot(dpS <- DivPart(q, MCS))
BetaS <- dpS$TotalBetaDiversity

# Temps
plot(MCT1 <- MetaCommunity(dfMC[, c(1,3)]))
plot(dpT1 <- DivPart(q, MCT1))
BetaT1 <- dpT1$GammaDiversity/dpT1$TotalAlphaDiversity
plot(MCT2 <- MetaCommunity(dfMC[, c(2,4)]))
plot(dpT2 <- DivPart(q, MCT2))
BetaT2 <- dpT2$GammaDiversity/dpT2$TotalAlphaDiversity
BetaTS <- sqrt(BetaT1*BetaT2)
plot(MergeMC(list(MCT1, MCT2)) -> MCT)
plot(dpT <- DivPart(q, MCT))
BetaT <- dpT$TotalBetaDiversity

# Décomposition
Beta
sqrt(BetaS)
sqrt(BetaST)
sqrt(BetaT)
sqrt(BetaTS)
sqrt(BetaS*BetaST*BetaT*BetaTS)

# BetaST original
Beta/BetaS/BetaT

# Diagramme euclidien
plot(c(0, sqrt(BetaT))~c(0, sqrt(BetaS)), type="n")
arrows(0, 0, sqrt(BetaS), sqrt(BetaT))
arrows(0, 0, sqrt(BetaS), 0)
arrows(0, 0, 0, sqrt(BetaT))

# Histogramme Rényi
barplot(cbind(c(log(BetaS), log(BetaT), log(BetaST), log(BetaTS)), NULL))
