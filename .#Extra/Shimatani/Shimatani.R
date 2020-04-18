library(ade4)         # pour phylog et taxo
library(cluster)      # pour l'arbre fonctionnel
library(entropart)    # pour les calculs directs de Simpson et Rao
library(mice)         # pour l'arbre fonctionnel
library(spatstat)     # pour les ppp et Kest
library(dbmss)        # pour les simulations


######################################################################################################################
## CONSTRUCTION DES ARBRES Taxonomy ET Functional
################################################

# Taxonomie
dfTaxo <- read.csv2(file="Taxonomie.csv", row.names=1)
Taxo <- as.taxo(dfTaxo)
TaxoPhy <- taxo2phylog(Taxo, add.tools=TRUE)
plot(TaxoPhy)
# Transformation de phylog en hclust
Taxonomy <- hclust(TaxoPhy$Wdist^2/2, "average")
plot(Taxonomy, h=-1)

# building individual trait database                                                            
trait1<- read.csv2(file="ind_traits_MICE.csv", header=T, sep=",", dec=".")  # leaf & stem economics
trait1<-trait1[,c(1:3,5,9,25)]                                              # sla & wsg
trait1$Name<- paste(trait1$Genus, trait1$species)                           # create vector to merge
trait2<- read.csv2(file="DataLifeTraits.csv", header=T)                     # life-history traits
trait2<-trait2[,c(3,6,10)]
trait<-merge(trait1, trait2, all.x=T, all.y=F)
trait$code<-paste(trait$Family, trait$Genus, trait$species, sep="_") 

# building species database    
Species <- data.frame(code=names(TaxoPhy$leaves))       
Species<-merge(Species, data.frame(code=names(tapply(trait$Hauteur, trait$code, median)),
                                   Height=as.numeric(tapply(trait$Hauteur, trait$code, median)),
                                   WSG=as.numeric(tapply(trait$sapwood_dens, trait$code, median)),
                                   SLA=as.numeric(tapply(trait$log_SLA, trait$code, median)),
                                   Seed=as.numeric(tapply(trait$Masse, trait$code, median))
                                   ), all.x=T, all.Y=F)
Species<-complete(mice(Species))                                            # using mice to fill gaps as in Paine-Baraloto
rownames(Species) <- Species$code
Species$code <- NULL
Species$Seed <- as.ordered(Species$Seed)                                    # seed mass, ordered factor

# functional tree using gower distance at the species level (all species including absent ones?)
Functional <- as.hclust(agnes(daisy(Species, metric="gower"), diss=TRUE, method="ward"))
plot(Functional)    


######################################################################################################################
# LECTURE DU JEU DE POINTS Paracou15 ET CONSTRUCTION D'UN JEU SIMPLIFIE X
###########################

ppTable<-read.csv2("Paracou 15.csv")
Paracou15 <- as.ppp(ppTable[ , c("X", "Y")], owin(c(min(ppTable[,"X"]), max(ppTable[,"X"])), c(min(ppTable[,"Y"]), max(ppTable[,"Y"]))))
Paracou15$marks <-  ppTable[ , "PointType"]
plot(Paracou15)

# Semis de points simplifié
X <- rmpoispp(lambda = c(.5, .3, .1, .05, .05)/1000, win= owin(c(0,250), c(0,250)), types=sample(names(TaxoPhy$leaves), 5))
plot(X)

######################################################################################################################
# SHIMATANI NEUTRE
###########################

Shimatani <- function(X, r, correction = "isotropic", Biased = TRUE) {
  
  fvCorrection <- function(x) {
    switch(correction,
           "isotropic" = x$iso,
           "translate" = x$trans,
           "none" = x$un
    )
  }
  
  # Summary
  Ns <- tapply(X$marks, X$marks, length)
  Ns <- Ns[!is.na(Ns)]
  Nall <- sum(Ns)
  Ps <- Ns/Nall
  
  # K all points
  Kall <- fvCorrection(Kest(X, r=r, correction=correction))
  
  # The point pattern is separated into a list of ppp for each mark
  pppList <- split(X, as.factor(X$marks))
  # K for each ppp
  KList <- lapply(pppList, Kest, r=r, correction=correction)
  Ks <- as.data.frame(lapply(KList, fvCorrection))
  # Ks is NA for species with a a single point. Should be 0
  Ks[is.na(Ks)] <- 0
  
  # Result
  Shr <- (1 - rowSums((Ks*rep(Ns*(Ns-1), each=dim(Ks)[1]))/(Kall*Nall*(Nall-1)))) * ifelse(Biased, 1, (Nall-1)/Nall)
  
  # Build a dataframe with r, theoretical value and S(r)
  ShiEstimate <- data.frame(r, Simpson(Ps), Shr)
  colnames(ShiEstimate) <- c("r", "Simpson", "S")
  
  # Return the values of Shimatani(r)
  return (fv(ShiEstimate, argu="r", ylab=quote(Shimatani(r)), valu="S", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "Simpson", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Simpson", "Estimated S(r)"), unitname=X$window$unit, fname="S"))

}

ShimataniEnvelope <-
function(X, r, NumberOfSimulations = 100, Alpha = 0.05, 
         correction = "isotropic", Biased = TRUE, SimulationType = "RandomLabeling", Global = FALSE) {
  
  # CheckdbmssArguments()
  
  # Choose the null hypothesis
  SimulatedPP <- switch (SimulationType,
                         RandomPosition = expression(rRandomPositionK(X, CheckArguments = FALSE)),
                         RandomLabeling = expression(rRandomLabeling(X, CheckArguments = FALSE))
  )
  if (is.null(SimulatedPP))
    stop(paste("The null hypothesis", sQuote(SimulationType), "has not been recognized."))
  # local envelope, keep extreme values for lo and hi (nrank=1)
  Envelope <- envelope(X, fun=Shimatani, nsim=NumberOfSimulations, nrank=1,
                       r=r, correction=correction, Biased=Biased,
                       #CheckArguments = FALSE,
                       simulate=SimulatedPP, savefuns=TRUE
  )
  attr(Envelope, "einfo")$H0 <- switch (SimulationType,
                                        RandomPosition = "Random Position",
                                        RandomLabeling = "Random Labeling"
  )
  # Calculate confidence intervals
  Envelope <- FillEnveloppe(Envelope, Alpha, Global)
  # Return the envelope
  return (Envelope)
}


# Applications
##############

# Semis de points sans structure
X <- rmpoispp(lambda = c(.5, .3, .1, .05, .05)/1000, win= owin(c(0,250), c(0,250)), types=sample(names(TaxoPhy$leaves), 5))
R=0:250
plot(Shimatani(X, R, "none", Biased=FALSE))
# Intervalle de confiance
ShEnv <- ShimataniEnvelope(X, R, NumberOfSimulations = 100, Alpha = 0.05, correction = "none", Biased = FALSE)
plot(ShEnv)

# Paracou 15
R <- c(seq(0, 10, 1), seq(15, 50, 5), seq(60, 120, 10))
plot(Shimatani(Paracou15, R, "none", Biased=FALSE))
# Intervalle de confiance
ShEnv <- ShimataniEnvelope(Paracou15, R, NumberOfSimulations = 10, Alpha = 0.05, correction = "none", Biased = FALSE)
plot(ShEnv)


######################################################################################################################
# SHIMATANI PHYLOGENETIQUE
###########################

# A distance r
PhyloShimatanir <- function(R, X, Tree = NULL, Normalized = FALSE, Biased = TRUE, ProgressBar = NULL) {
  
  # Transform the tree into a distance matrix
  if (is.null(Tree)){
    Dist <- 1-diag(length(levels(X$marks)))
    rownames(Dist) <- colnames(Dist) <- levels(X$marks)
    Tree <- hclust(dist(Dist/sqrt(2)))
  } else {
    Dist <- as.matrix(cophenetic(Tree))
    if (Normalized)
      Dist <- Dist/max(Dist)
  }
  
  # In each point neighborhood, count pairs and sum the distances
  nbdPairs <- function(Y, current, ...) {
    ReferenceMark <- as.character(current$marks)
    Distances <- lapply(as.character(Y$marks), function(marks) Dist[ReferenceMark, marks])
    c(length(Distances), sum(unlist(Distances)))
  }
  if (!is.null(ProgressBar))
    setTxtProgressBar(ProgressBar, R)
  nbdSummary <- rowSums(applynbd(X, nbdPairs, R=R, exclude=Biased))
  nbdSummary[2]/nbdSummary[1]
}

# Pour un vecteur de distances, avec barre de progression
PhyloShimatani <- function(X, r, Tree = NULL, Normalized = FALSE, Biased = TRUE, CalculateRao = TRUE, ShowProgressBar = TRUE) {
  # ProgressBar
  if (ShowProgressBar) {
    ProgressBar <- txtProgressBar(min=min(r), max=max(r))
  } else {
    ProgressBar <- NULL
  }
  # Calculate pS(r)
  pSh <- sapply(r, PhyloShimatanir, X, Tree, Normalized, Biased, ProgressBar)
  # ProgressBar
  if (!is.null(ProgressBar))
    close(ProgressBar)
  # Build a dataframe with r, theoretical value and S(r)
  Ns <- tapply(X$marks, X$marks, length)
  Nall <- sum(Ns)
  Ps <- Ns/Nall
  ShiEstimate <- data.frame(r, ifelse(CalculateRao, Rao(Ps, Tree), NA), pSh)
  colnames(ShiEstimate) <- c("r", "Rao", "pS")
  
  # Return the values of Shimatani(r)
  return (fv(ShiEstimate, argu="r", ylab=quote(PhyloShimatani(r)), valu="pS", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "Rao", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Rao", "Estimated pS(r)"), unitname=X$window$unit, fname="pS"))
}


# Applications
##############

# Semis de points sans structure
X <- rmpoispp(lambda = c(.5, .3, .1, .05, .05)/1000, win= owin(c(0,250), c(0,250)), types=sample(names(TaxoPhy$leaves), 5))
R=seq(0, 100, 5)
pShi <- PhyloShimatani(X, R, Taxonomy)
plot(pShi)
# Calcul de Rao avec divc
divc(as.data.frame(tapply(X$marks, X$marks, length)), as.dist(as.matrix(TaxoPhy$Wdist)[levels(X$marks),levels(X$marks)]))
# Intervalle de confiance
pShEnv <- PhyloShimataniEnvelope(X, R, NumberOfSimulations = 100, Alpha = 0.05, Tree=Taxonomy)
plot(pShEnv)


# Paracou 15, taxonomie (long !)
R <- c(seq(0, 10, 1), seq(15, 50, 5), seq(60, 120, 10))
pShi <- PhyloShimatani(Paracou15, R, Taxonomy)
plot(pShi)
# Intervalle de confiance (très long !)
pShEnv <- PhyloShimataniEnvelope(X, R, NumberOfSimulations = 2, Alpha = 0.05, Tree=Taxonomy)
plot(pShEnv)
# Paracou 15, fonctionnel pour comparaison (long !)
pShi <- PhyloShimatani(Paracou15, R, Functional)
plot(pShi)

######################################################################################################################
# SHIMATANI PHYLO avec correction des effets de bord
###########################

PhyloShimataniK <- function(X, r, Tree = NULL, Normalized = FALSE, correction = "isotropic", Biased = TRUE, CalculateRao = TRUE, ShowProgressBar = TRUE) {
  
  fvCorrection <- function(x) {
    switch(correction,
           "isotropic" = x$iso,
           "translate" = x$trans,
           "none" = x$un
    )
  }
  
  # Transform the tree into a distance matrix
  if (is.null(Tree)){
    Dist <- 1-diag(length(levels(X$marks)))
    rownames(Dist) <- colnames(Dist) <- levels(X$marks)
    Tree <- hclust(dist(Dist/sqrt(2)))
  } else {
    Dist <- as.matrix(cophenetic(Tree))
    if (Normalized)
      Dist <- Dist/max(Dist)
  }
  
  # Summary
  Ns <- tapply(X$marks, X$marks, length)
  Nall <- sum(Ns)
  Ps <- Ns/Nall
  
  # K all points
  Kall <- fvCorrection(Kest(X, r=r, correction=correction))
  
  # Calculate all bivariate K's
  PointTypes <- levels(marks(X))
  # ProgressBar
  TxtProgressBarValue <- 0
  if (ShowProgressBar) {
    ProgressBar <- txtProgressBar(min=0, length(PointTypes)^2)
  } else {
    ProgressBar <- NULL
  }
  
  WeightedSumOfDistances <- vector(mode="numeric", length=length(r))
  for (i in PointTypes) {
    for (j in PointTypes) {
      Kij <- Kcross(X, i, j, r=r, correction=correction)
      KijVector <- fvCorrection(Kij)
      # K is NA for species with a a single point. Should be 0
      KijVector[is.na(KijVector)] <- 0
      WeightedSumOfDistances <- WeightedSumOfDistances + KijVector*Ns[i]*(Ns[j]-ifelse(i==j, 1, 0))*Dist[i,j]
      if (!is.null(ProgressBar)) {
        TxtProgressBarValue <- TxtProgressBarValue+1
        setTxtProgressBar(ProgressBar, TxtProgressBarValue)
      }
    }
  }
  
  # Result
  pSh <- WeightedSumOfDistances/(Kall*Nall*(Nall-1)) * ifelse(Biased, 1, (Nall-1)/Nall)
  ShiEstimate <- data.frame(r, ifelse(CalculateRao, AllenH(Ps, 2, Tree, Normalize=Normalized, CheckArguments=FALSE), NA), pSh)
  colnames(ShiEstimate) <- c("r", "Rao", "pS")

  # ProgressBar
  if (!is.null(ProgressBar))
    close(ProgressBar)
  
  # Return the values of Shimatani(r)
  return (fv(ShiEstimate, argu="r", ylab=quote(PhyloShimatani(r)), valu="pS", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "Rao", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Rao", "Estimated pS(r)"), unitname=X$window$unit, fname="pS"))

}

PhyloShimataniEnvelope <-
  function(X, r, NumberOfSimulations = 100, Alpha = 0.05, 
           Tree = NULL, Normalized = FALSE, correction = "isotropic", Biased = TRUE, SimulationType = "RandomLabeling", Global = FALSE, ShowProgressBar = FALSE) {
    
    #CheckdbmssArguments()
    
    # Choose the null hypothesis
    SimulatedPP <- switch (SimulationType,
                           RandomPosition = expression(rRandomPositionK(X, CheckArguments = FALSE)),
                           RandomLabeling = expression(rRandomLabeling(X, CheckArguments = FALSE))
    )
    if (is.null(SimulatedPP))
      stop(paste("The null hypothesis", sQuote(SimulationType), "has not been recognized."))
    # local envelope, keep extreme values for lo and hi (nrank=1)
    Envelope <- envelope(X, fun=PhyloShimataniK, nsim=NumberOfSimulations, nrank=1,
                         r=r, Tree=Tree, Normalized=Normalized, correction=correction, Biased=Biased,
                         CalculateRao=FALSE, ShowProgressBar=ShowProgressBar,
                         # CheckArguments = FALSE,
                         simulate=SimulatedPP, savefuns=TRUE
    )
    attr(Envelope, "einfo")$H0 <- switch (SimulationType,
                                          RandomPosition = "Random Position",
                                          RandomLabeling = "Random Labeling"
    )
    # Calculate confidence intervals
    Envelope <- FillEnveloppe(Envelope, Alpha, Global)
    # Return the envelope
    return (Envelope)
  }

######################################################################################################################
# VERIFICATION
###########################
# Semis de points sans structure
X <- rmpoispp(lambda = c(.5, .3, .1, .05, .05)/1000, win= owin(c(0,250), c(0,250)), types=sample(names(TaxoPhy$leaves), 5))
R=seq(0, 100, 5)
plot(X)

# PhiloShimatani calculé par applynbd = PhiloShimatani calculé par K sans correction des effets de bord
system.time(pShi <- PhyloShimatani(X, R, Taxonomy))
plot(pShi)
system.time(pShiK <- PhyloShimataniK(X, R, Taxonomy, correction="none"))
plot(pShiK)

# PhiloShimatani calculé par K avec un arbre dont toutes les distances valent 1 = Shimatani neutre
# Correction des effets de bord par défaut (isotropic)
pShiK0 <- PhyloShimataniK(X, R)
plot(pShiK0)
pShiNeutre <- Shimatani(X, R)
plot(pShiNeutre)

# Envelope
pShEnv <- PhyloShimataniEnvelope(X, R, NumberOfSimulations=100, Alpha=0.05, Tree=Taxonomy)
plot(pShEnv)
# Rao
abline(h=Rao(tapply(X$marks, X$marks, length)/X$n,Taxonomy), lty=4, col="blue")
