# Rao Entropy along a vector of distances
RaoR <- function(X, r, Tree = NULL, Normalized = FALSE, correction = "isotropic", 
              Biased = TRUE, CalculateRao = TRUE, ShowProgressBar = TRUE, CheckArguments = TRUE)
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
  Ns <- Ns[!is.na(Ns)]
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
