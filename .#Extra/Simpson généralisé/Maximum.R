##### Recherche de la distribution de probabilité qui maximise zeta
GenSimpsonOpp <- function(Ps) -GenSimpson(Ps, r, FALSE)
Constraint <- function(Ps) sum(Ps)
library(Rsolnp) 

# Paramètres
r <- 100
S <- 97

# Valeurs initiales pour la résolution de l'équation
(Ps <- c(.5, rep(.5/(S-1), S-1)))
# Résolution
Solution <- solnp(Ps, fun=GenSimpsonOpp, eqfun=Constraint, eqB=1, LB=rep(0, S), UB=rep(1, S))

# Solution et valeur de la fonction
Solution$pars
GenSimpson(Solution$pars, r)
GenSimpsonD(Solution$pars, r)

# Solution si r >> S
1-(S-1)/(r+1)
1/(r+1)
# Valeur de Simpson généralisé
GenSimpson(c(1-(S-1)/(r+1), rep(1/(r+1), S-1)), r)
GenSimpsonD(c(1-(S-1)/(r+1), rep(1/(r+1), S-1)), r)

# Optimisation de f (Simpson généralisé avec toutes probas égales sauf 1)
Simpsonrx <- function(x) x*(1-x)^r +(1-x)*(1-(1-x)/(S-1))^r
# Valeur de la probabilité libre qui maximise f
optim(1/S, Simpsonrx, method = "Brent", lower=0, upper=1, control=list(fnscale=-1))$par
