# Données : distribution lognormale
Nspecies <- 100
sdlog <- 2
Size <- 400

rCommunity(n = 1, size = Size, S=Nspecies, Distribution = "lnorm", sd=sdlog) -> NsRef
plot(NsRef, Distribution="lnorm")


### Profil de zeta
r.seq <- 1:(Size-1)
cpzeta <- CommunityProfile(FUN = GenSimpson, NorP = as.ProbaVector(NsRef), q.seq=r.seq)
plot(cpzeta)

# Courbe d'accumulation
plot(c(0, r.seq), c(0, cumsum(cpzeta$y)), type="l", xlab="Sampling Effort", ylab="Expected Number of Species")

# Nombre effectif
cpzetaD <- CommunityProfile(FUN = GenSimpsonD, NorP = as.ProbaVector(NsRef), q.seq=r.seq)
plot(cpzetaD)
abline(h=Nspecies)
# Quelle valeur de r donne le nombre d'espèces ?
cpzetaD$x[which.max(cpzetaD$y > Nspecies)]
abline(v=cpzetaD$x[which.max(cpzetaD$y > Nspecies)])

### Estimateur z
library("EntropyEstimation")
N <- sum(NsRef)
y <- as.numeric(sapply(r.seq, function(r) GenSimp.z(NsRef, r)))
ysd <- as.numeric(sapply(r.seq, function(r) GenSimp.z(NsRef, r)))
Zenvelope <- as.CommunityProfile(r.seq, y, y-2*ysd/sqrt(N), y+2*ysd/sqrt(N))
plot(Zenvelope)
CEnvelope(cpzeta, col="Green")

#### Estimation du nombre d'espèces: équivalent à la correction de ZG
NumberOfSimulations <- 5
Size <- 10000
r.seq <- 1:(Size-1)
MCSim <- rCommunity(NumberOfSimulations, size=Size, NorP=NsRef)
zetaP <- apply(MCSim$Nsi, 2, function(Nsi) as.numeric(sapply(r.seq, function(r) GenSimp.z(Nsi, r))))
(EstS <- colSums(zetaP))
mean(EstS)

#### Comparaison avec HCDT
cpHCDT <- CommunityProfile(FUN = Diversity, NorP = as.ProbaVector(NsRef), q.seq=seq(0, 5, .1))
plot(cpHCDT)

#### Comparaison avec Hurlbert
cpHurlbert <- CommunityProfile(FUN = HurlbertD, NorP = NsRef, q.seq=r.seq[-1])
plot(cpHurlbert)
