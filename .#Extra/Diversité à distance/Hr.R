# Entropy along a vector of distances
Hr <- function(X, r, nPoints, CenterOnPoints = FALSE, FUN, ..., Limit = TRUE, ShowProgressBar = TRUE, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()  

  if (CenterOnPoints) {
    # Keep nPpoints of the ppp for centers
    Centers <- rthin(X, nPoints/X$n)
  } else {
    # Draw centers along a regular grid
    Centers <- rsyst(win=X$window, dx=sqrt(area.owin(X)/nPoints))  
  }
  
  # ProgressBar
  if (ShowProgressBar) {
    ProgressBar <- txtProgressBar(min=0, max=length(r)*Centers$n)
    pbValue <- 0
  } else {
    ProgressBar <- NULL
  }

  Alphar <- vector(mode="numeric", length=length(r))
  Alphai <- ni <- vector(mode="numeric", length=Centers$n)
  for (j in r) {
    for (i in 1:Centers$n) {
      # Calculate the function at each center point
      de <- DiscEntropy(j, c(Centers$x[i], Centers$y[i]), X, FUN, ...)
      Alphai[i] <- de$H
      ni[i] <- de$n-1
      if (!is.null(ProgressBar))
        setTxtProgressBar(ProgressBar, pbValue)
        pbValue <- pbValue+1
    }
    # Average it to get the value at distance r
# print(Alphai)
    Alphar[which(r==j)] <- sum(Alphai*ni, na.rm=TRUE)/sum(ni[!is.na(Alphai)])
  
  }
  # ProgressBar
  if (!is.null(ProgressBar))
    close(ProgressBar)

  # Limit value
  if (Limit) {
    # Calculate the entropy of a disc containing the whole window
    LimitValue <- DiscEntropy(diameter(X$window), c(Centers$x[1], Centers$y[1]), X, FUN, ...)
  }
  # Build a dataframe with r and Hr(r)
  ShiEstimate <- data.frame(r, ifelse(Limit, LimitValue, NA), Alphar)
  colnames(ShiEstimate) <- c("r", "Limit", "Hr")
  
  # Return the values of Shimatani(r)
  return (fv(ShiEstimate, argu="r", ylab=quote(Hr(r)), valu="Hr", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "Limit", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Limit Value", "Estimated Hr(r)"), unitname=X$window$unit, fname="Hr"))
}
