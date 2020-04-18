library(entropart)    # pour les calculs directs de Simpson et Rao
library(dbmss)        # pour les simulations


######################################################################################################################
# LECTURE DU JEU DE POINTS Paracou15 ET CONSTRUCTION D'UN JEU SIMPLIFIE X
###########################

ppTable<-read.csv2("Paracou 15.csv")
Paracou15 <- as.ppp(ppTable[ , c("X", "Y")], owin(c(min(ppTable[,"X"]), max(ppTable[,"X"])), c(min(ppTable[,"Y"]), max(ppTable[,"Y"]))))
Paracou15$marks <- ppTable[ , "PointType"]
plot(Paracou15)
# Taxonomie
dfTaxo <- read.csv2(file="Taxonomie.csv", row.names=1)
Taxo <- as.taxo(dfTaxo)
TaxoPhy <- taxo2phylog(Taxo, add.tools=TRUE)
plot(TaxoPhy)

# un carré
P15 <- Paracou15[owin(c(0,125), c(0,125))]
# Semis de points simplifié
X <- rmpoispp(lambda = c(.5, .3, .1, .05, .05)/1000, win= owin(c(0,250), c(0,250)), types=sample(names(TaxoPhy$leaves), 5))
plot(X)

#########################################################
# Carré 4 points
win <- owin(c(0,10), c(0,10))
df <- data.frame(x=c(1,2.1,1, 8,9.1, 2), y=c(1,2.1,2.1, 8,9.1, 1), PointType=c("A", "A", "B", "A", "B", "C"))
X4 <- as.wmppp(df, win)
plot(X4)
q=2
R = seq(1,3,1)
FUN=bcSimpson
Hr(X4, R, nPoints=X4$n, CenterOnPoints=TRUE, FUN=FUN, Correction="Lande") -> Entropyr
plot(Entropyr)
SR <- SimpsonR(X4, c(0, R), correction="none")
plot(SR)
#########################################################
# Carré 50 points
win <- owin(c(0,10), c(0,10))
X <- rmpoint(n=10, types=c("A", "B", "C", "D"), win=win)
plot(X)
q=2
R = seq(1,5,1)
FUN=bcSimpson
Hr(X, R, nPoints=X$n, CenterOnPoints=TRUE, FUN=FUN, Correction="Lande") -> Entropyr
plot(Entropyr)
Entropyr$Hr
SR <- SimpsonR(X, c(0, R), correction="none")
plot(SR)
SR$Sr

#########################################################
# Calcul

q=2
R=seq(5, 50, 5)
# Tsallis
FUN=bcTsallis
Hr(P15, R, nPoints=P15$n, CenterOnPoints=TRUE, FUN=FUN, q=q, Correction="None") -> Entropyr
plot(Entropyr)
# Shimatani
SR <- SimpsonR(P15, c(0, R), correction = "none", Biased = F)
plot(SR)

HrEnvelope(Paracou15, R, nPoints=1, FUN=bcTsallis, q=q,  Correction="None", NumberOfSimulations = 4) -> Envelopper
plot(Envelopper)
# Diversité
Order <- q
plot(Envelopper, expq(., Order) ~ r)