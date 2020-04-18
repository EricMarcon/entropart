# Simpson Entropy along a vector of distances
SimpsonR <- function(X, r, correction = "isotropic", Biased = FALSE, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()

  # Get the appropirate value of the fv returned by spatstat according to the correction
  fvCorrection <- function(x) {
    switch(correction,
           "isotropic" = x$iso,
           "translate" = x$trans,
           "none" = x$un
    )
  }
  
  # Count Ns
  if (is.wmppp(X)) {
    Species <- marks(X)$PointType
  } else {
    Species <- marks(X)
  }
  Ns <- tapply(Species, Species, length)
  Ns <- Ns[!is.na(Ns)]
  Nall <- sum(Ns)
  Ps <- Ns/Nall
  
  # K all points
  Kall <- fvCorrection(Kest(X, r=r, correction=correction))
  
  # The point pattern is separated into a list of ppp for each mark
  if (is.wmppp(X)) {
    pppList <- split(X, as.factor(X$marks$PointType))
  } else {
    pppList <- split(X, as.factor(X$marks))
  }
  # K for each ppp
  KList <- lapply(pppList, Kest, r=r, correction=correction)
  Ks <- as.data.frame(lapply(KList, fvCorrection))
  # Ks is NA for species with a a single point. Should be 0
  Ks[is.na(Ks)] <- 0
  
  # Result
  Shr <- (1 - rowSums((Ks*rep(Ns*(Ns-1), each=dim(Ks)[1]))/(Kall*Nall*(Nall-1)))) * ifelse(Biased, 1, (Nall-1)/Nall)
  
  # Build a dataframe with r, theoretical value and S(r)
  ShiEstimate <- data.frame(r, ifelse(Biased, Simpson(Ps), bcSimpson(Ns)), Shr)
  colnames(ShiEstimate) <- c("r", "Simpson", "Sr")
  
  # Return the values of SimpsonR(r)
  return (fv(ShiEstimate, argu="r", ylab=quote(Sr(r)), valu="Sr", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "Simpson", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Simpson", "Estimated Sr(r)"), unitname=X$window$unit, fname="Sr"))
  
}
