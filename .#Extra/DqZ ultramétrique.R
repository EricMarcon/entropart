library("entropart")
############### Données ####
# q
q <- 2
# Matrice de distances, 3 espèces
D <- matrix(c(  0, 2, 3,
                2, 0, 3,
                3, 3, 0), nrow=3, byrow = TRUE)
colnames(D) <- rownames(D) <- c("A", "B", "C")
Tree <- hclust(as.dist(D))
Z <- 1 - D/max(D)

# Fréquences: 2 communautés
Nsi <- matrix(c(5, 3,
                5, 8,
                5, 1), nrow=3, byrow = TRUE)
rownames(Nsi) <- c("A", "B", "C")
# Poids
Wi <- c(2/3, 1/3)

# Création de la MC
MC <- MetaCommunity(Nsi, Weights=Wi)

# Diversité de L&C
Dqz(MC$Ps, q, Z)

# Différent, si q<>2, de
# Approche alternative (Proposition A7)
Zt <- matrix(c( 1, 1, 0, 0, 0, 
                1, 1, 1, 1, 0,
                0, 0, 1, 1, 0,
                1, 1, 1, 1, 0,
                0, 0, 0, 0, 1), nrow=5, byrow = TRUE)

Pt <- c(2/3*MC$Ps["A"], 1/3*MC$Ps["A"], 2/3*MC$Ps["B"], 1/3*MC$Ps["B"], 1*MC$Ps["C"])
Dqz(Pt, q, Zt)


# Identique à :
PhyloDiversity(MC$Ps, q, Tree)$Total

