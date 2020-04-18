##### PCoA d'ordre q ###
q <- 1
# Deux communautés de 3 espèces
mcmat <-matrix(as.integer(1+runif(6)*10), ncol=2) # Aléatoire sans 0
plot(MC <- MetaCommunity(mcmat))
plot(DP <- DivPart(q, MC))

# Distances entre espèces et c1
Dist <- matrix(1, nrow=4, ncol=4)
diag(Dist) <- 0

# Distances
Dist[1,4] <- Dist[4,1] <- lnq(1/MC$Psi[1], q)
Dist[2,4] <- Dist[4,2] <- lnq(1/MC$Psi[2], q)
Dist[3,4] <- Dist[4,3] <- lnq(1/MC$Psi[3], q)

library("ade4")
dDist <- as.dist(sqrt(Dist))
is.euclid(dDist)

edDist <- lingoes(dDist)
scatter(dudi.pco(edDist, scannf = F))
