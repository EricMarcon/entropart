# Valeurs des probabilités. Pour zoomer sur le pic du fond, utiliser seq(0, .1, .001)
ps <- seq(0, 1, .01)
# q : paramètre, à modifier
q <- 4

# S, fixe
S <- 3

# Opposée (l'optimisation ne sait que minimiser) de la fonction (Simpson généralisé à l'ordre q)
Simpsonq <- function(Ps) -sum(Ps*(1-Ps)^q)

# Utilitaire : renvoie un vecteur des valeurs de la fonction pour une valeur de p choisie, 
# toutes les valeurs possibles de la deuxième probabilité, la troisième est la différence
# Nécessaire pour le graphique en 3D
Simpsonp1 <- function(p1, q) {
  Result <- sapply(ps, function(p2) p1*(1-p1)^q + p2*(1-p2)^q + (1-p1-p2)*(p1+p2)^q)
  # Elimination des combinaisons de probabilités dont la somme est > 1
  Result[1-p1-ps<=0] <- NA 
  return(Result)
}

# Valeurs de la fonction pour toutes les valeurs de p. Renvoie une matrice de valeurs pour p x p
f <- sapply(ps, Simpsonp1, q)

# Graphique 3D
library(plot3D)
par(mar = c(0, 2, 0, 2))
persp3D(z=f, x=ps, y=ps, facets=FALSE, curtain=TRUE, xlab="p1", ylab="p2", zlab="Simpson", ticktype="detailed", colkey=list(length=.5), theta = 135) -> Projection
# Points : triangle pour toutes les p égales
points(trans3d(1/S, 1/S, (1-1/S)^q, pmat=Projection), col="black", pch=17)
# 3 points pour deux proba égales à 1/(q+1) et une égale à 1-S/(q+1) : diamants
points(trans3d(1/(q+1), 1/(q+1), -Simpsonq(c(rep(1/(q+1), 2), 1-2/(q+1))), pmat=Projection), col="black", pch=5)
points(trans3d(1/(q+1), 1-2/(q+1), -Simpsonq(c(rep(1/(q+1), 2), 1-2/(q+1))), pmat=Projection), col="black", pch=5)
points(trans3d(1-2/(q+1), 1/(q+1), -Simpsonq(c(rep(1/(q+1), 2), 1-2/(q+1))), pmat=Projection), col="black", pch=5)

# Recherche du maximum.
# Valeur de départ, une proba à 1/2 et les autres égales
(Ps <- c(.5, rep(.5/(S-1), S-1)))
# Contrainte : somme des p = 1
Constraint <- function(Ps) sum(Ps)
# Résolution
library(Rsolnp)
(Solution <- solnp(Ps, fun=Simpsonq, eqfun=Constraint, eqB=1, LB=rep(0, S), UB=rep(1, S)))
# Valeurs des proba optimales
Solution$pars
# Tracé d'un point optimal (croix) sur le pic du fond
points(trans3d(Solution$pars[2], Solution$pars[3], -Simpsonq(Solution$pars), pmat=Projection), col="black", pch=3)

