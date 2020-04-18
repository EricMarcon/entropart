library(entropart)
# 3 communautes log-normales
C1 <- rCommunity(1, size=300)
C2 <- rCommunity(1, size=300)
C3 <- rCommunity(1, size=300)
# Poids arbitraires
w <- 1:3 # rep(1, 3)
# Normalisation des poids
w <- w/sum(w)
# Création d'une métacommunauté
MC <-   MetaCommunity(cbind(C1, C2, C3), Weights = w)
# Métacommunautés réduites à deux communautés
C12 <-  MetaCommunity(cbind(C1, C2), Weights = w[-3])
C23 <-  MetaCommunity(cbind(C2, C3), Weights = w[-1])
C13 <-  MetaCommunity(cbind(C1, C3), Weights = w[-2])
# Métacommunautés rassemblant les précédentes et la communauté complémentaire
MC_1 <- MetaCommunity(cbind(C23$Ns, C1), Weights = c(sum(w[-1]), w[1]))
MC_2 <- MetaCommunity(cbind(C13$Ns, C2), Weights = c(sum(w[-2]), w[2]))
MC_3 <- MetaCommunity(cbind(C12$Ns, C3), Weights = c(sum(w[-3]), w[3]))
# Ordre de la diversité
q <- 0.5
# Diversité beta des métacommunautés
DPmc <- DivPart(q, MC)
DP12 <- DivPart(q, C12)
DP23 <- DivPart(q, C23)
DP13 <- DivPart(q, C13)
DP_1 <- DivPart(q, MC_1)
DP_2 <- DivPart(q, MC_2)
DP_3 <- DivPart(q, MC_3)

# Vérification de la décomposition: somme égale à 3 diversités beta
DP12$TotalBetaEntropy*sum(w[-3]) + DP23$TotalBetaEntropy*sum(w[-1]) + DP13$TotalBetaEntropy*sum(w[-2]) +
  DP_1$TotalBetaEntropy + DP_2$TotalBetaEntropy + DP_3$TotalBetaEntropy
3*DPmc$TotalBetaEntropy

# Diversité beta_i entre l'assemblage de deux communautés et la métacommunauté
beta1_1 <- (DP_1$GammaEntropy - DP_1$CommunityAlphaEntropies)[1]
beta2_1 <- (DP_2$GammaEntropy - DP_2$CommunityAlphaEntropies)[1]
beta3_1 <- (DP_3$GammaEntropy - DP_3$CommunityAlphaEntropies)[1]
# Diversité beta_i entre la communauté non regroupée et la métacommunauté
beta1_2 <- (DP_1$GammaEntropy - DP_1$CommunityAlphaEntropies)[2]
beta2_2 <- (DP_2$GammaEntropy - DP_2$CommunityAlphaEntropies)[2]
beta3_2 <- (DP_3$GammaEntropy - DP_3$CommunityAlphaEntropies)[2]


# Vérification de la décomposition: somme égale à 2 diversités beta
MC_1$Wi[1]*(beta1_1 + DP23$TotalBetaEntropy) + 
  MC_2$Wi[1]*(beta2_1 + DP13$TotalBetaEntropy) + 
  MC_3$Wi[1]*(beta3_1 + DP12$TotalBetaEntropy)
2*DPmc$TotalBetaEntropy

# Métacommunauté des assemblages deux à deux.
MC6C <- MetaCommunity(cbind(C12$Ns, C13$Ns, C23$Ns), Weights = c(sum(w[-3]), sum(w[-2]), sum(w[-1])))
# Vérification de la décomposition
MC_1$Wi[1]*beta1_1 + MC_2$Wi[1]*beta2_1 + MC_3$Wi[1]*beta3_1 
2*DivPart(q, MC6C)$TotalBetaEntropy

# Matrice de distance deux à deux
D <- matrix(c(0,                     DP12$TotalBetaEntropy, DP13$TotalBetaEntropy,
              DP12$TotalBetaEntropy, 0,                     DP23$TotalBetaEntropy, 
              DP13$TotalBetaEntropy, DP23$TotalBetaEntropy, 0 ), nrow = 3)
library(ade4)
is.euclid(as.dist(D))
scatter(dudi.pco(as.dist(D), scannf = F))

# Ajout de la MC
Dmc <- matrix(c(0,                       2*DP12$TotalBetaEntropy, 2*DP13$TotalBetaEntropy, beta1_2, 
                2*DP12$TotalBetaEntropy, 0,                       2*DP23$TotalBetaEntropy, beta2_2,
                2*DP13$TotalBetaEntropy, 2*DP23$TotalBetaEntropy, 0,                       beta3_2,
                beta1_2,                 beta2_2,                 beta3_2,                 0 ), nrow = 4)
is.euclid(as.dist(Dmc))
scatter(dudi.pco(as.dist(Dmc), scannf = F), posieig = "none")

# Distance au centre d'une paire de communautés : bof
MC1_12 <- MetaCommunity(cbind(C1, C12$Ns), Weights = rep(1, 2))
DP1_12 <- DivPart(q, MC1_12)
MC2_12 <- MetaCommunity(cbind(C2, C12$Ns), Weights = rep(1, 2))
DP2_12 <- DivPart(q, MC2_12)
D12 <- matrix(c(0,                     DP12$TotalBetaEntropy, DP1_12$TotalBetaEntropy,
                DP12$TotalBetaEntropy, 0,                     DP2_12$TotalBetaEntropy, 
                DP1_12$TotalBetaEntropy, DP2_12$TotalBetaEntropy, 0 ), nrow = 3)
is.euclid(as.dist(D12))


# Distance au centre d'une paire de communautés : décomposition

