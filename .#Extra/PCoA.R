library("entropart")

library("vegan")
data(BCI)

# 2 communautés de 3 espèces
BCI[2,50] <- 2
BCI[2,49] <- 20
mcBCI <- MetaCommunity(t(BCI[1:3, 49:51]))
plot(mcBCI)

# Affichage dans le plan des deux espèces
Psi <- mcBCI$Psi
RDA <- rda(Psi)
RDA$CA$rank
RDA$CA$eig
plot(RDA)

# Diversité beta
q <- 2
# Entropie  des communautés
(gamma <- Tsallis(mcBCI$Ps, q))
(beta <- apply(mcBCI$Psi, 2, function(p) TsallisBeta(p, mcBCI$Ps, q)))
(alpha <- apply(mcBCI$Psi, 2, function(p) Tsallis(p, q)))
sum((alpha+beta) %*% mcBCI$Wi)

# Distance deux à deux
BetaDist <- function(Ps1, Ps2, W, q) {
  Ps <- matrix(c(Ps1, Ps2), ncol = 2) %*% (W/sum(W))
  crossprod(c(TsallisBeta(Ps1, Ps, q), TsallisBeta(Ps2, Ps, q)), W/sum(W))
}

# La divergence totale égale la div entre 1 et 2 + celle du groupe 12 + celle de 3
Ps12 <- matrix(c(mcBCI$Psi[,1], mcBCI$Psi[,2]), ncol = 2) %*% (mcBCI$Wi[c(1,2)]/sum(mcBCI$Wi[c(1,2)]))
rownames(Ps12) <- rownames(mcBCI$Psi)
BetaDist(mcBCI$Psi[, 1], mcBCI$Psi[, 2], mcBCI$Wi[c(1, 2)], q)*sum(mcBCI$Wi[c(1,2)]) 
  + BetaDist(Ps12, mcBCI$Psi[, 3], c(2/3, 1/3), q)
crossprod(beta, mcBCI$Wi)

# Divergences 2 à 2
BetaDist(mcBCI$Psi[, 1], mcBCI$Psi[, 2], mcBCI$Wi[c(1, 2)], q)
BetaDist(mcBCI$Psi[, 1], mcBCI$Psi[, 3], mcBCI$Wi[c(1, 3)], q)
BetaDist(mcBCI$Psi[, 2], mcBCI$Psi[, 3], mcBCI$Wi[c(2, 3)], q)


# Matrice de distance entre tous les points  : 3 espèces, 2 communautés, métacommunauté
# Réduction à deux communautés pour commencer
mcBCI <- MetaCommunity(t(BCI[1:2, 49:51]))
D <- matrix(rep(1, 36), ncol=6)
diag(D) <- 0
colnames(D) <- rownames(D) <- c("sp1", "sp2", "sp3", "C1", "C2", "MC")

# Distances entre espèces et communautés
D[1, 4] <- D[4, 1] <- TsallisBeta(c(1, 0, 0), mcBCI$Psi[, 1], q)
D[2, 4] <- D[4, 2] <- TsallisBeta(c(0, 1, 0), mcBCI$Psi[, 1], q)
D[3, 4] <- D[4, 3] <- TsallisBeta(c(0, 0, 1), mcBCI$Psi[, 1], q)
D[1, 5] <- D[5, 1] <- TsallisBeta(c(1, 0, 0), mcBCI$Psi[, 2], q)
D[2, 5] <- D[5, 2] <- TsallisBeta(c(0, 1, 0), mcBCI$Psi[, 2], q)
D[3, 5] <- D[5, 3] <- TsallisBeta(c(0, 0, 1), mcBCI$Psi[, 2], q)

# Distances entre espèces et méta-communauté
D[1, 6] <- D[6, 1] <- TsallisBeta(c(1, 0, 0), mcBCI$Ps, q)
D[2, 6] <- D[6, 2] <- TsallisBeta(c(0, 1, 0), mcBCI$Ps, q)
D[3, 6] <- D[6, 3] <- TsallisBeta(c(0, 0, 1), mcBCI$Ps, q)


# Distances entre communautés et méta-communauté
dp <- DivPart(q, mcBCI, Biased = T)
D[4, 6] <- D[6, 4] <- dp$CommunityBetaEntropies[1]
D[5, 6] <- D[6, 5] <- dp$CommunityBetaEntropies[2]

# Distance entre C1 et C2 
D[4, 5] <- D[5, 4] <- BetaDist(mcBCI$Psi[, 1], mcBCI$Psi[, 2], mcBCI$Wi[c(1, 2)], q)

Dc=D[-c(4,5), -c(4,5)]

# Matrice de distance: tests
library("ade4")
dD <- as.dist(sqrt(2*Dc))
is.euclid(dD)
scatter(dudi.pco(dD, scannf = F,  nf = 2))

# Placer C1 dans l'espace des espèces
library("ade4")
DC1 <- sqrt(as.dist(D[1:4, 1:4]))
is.euclid(DC1)
pcoa1 <- dudi.pco(DC1, scannf = F,  nf = 2)
scatter(pcoa1)

DC2 <- sqrt(as.dist(D[c(1:3, 5), c(1:3, 5)]))
is.euclid(DC2)
pcoa2 <- dudi.pco(DC2, scannf = F, nf = 2)
scatter(pcoa2)
