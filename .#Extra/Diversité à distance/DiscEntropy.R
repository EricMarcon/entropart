# Entropy in a circular point pattern
DiscEntropy <- function(r, xy, X, FUN, ..., CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()  

#Xb=X
#X=X[!((X$x == xy[1]) & (X$y == xy[2]))]
  # Create a disc ppp
  Disc <- spatstat:::disc(radius=r, centre=xy)
  # Intersect it with X
  Xr <- X[Disc]
#X=Xb
  # Eliminate point weights id X is a wmppp (defined in package dbmss)
  if (is.wmppp(X))
    marks(Xr) <- marks(Xr)$PointType
  # Count species anundances
  Ns <- tapply(Xr$marks, Xr$marks, length)
  # Eliminate NA due to species with abundance 0
  Ns <- Ns[!is.na(Ns)]
  if (sum(Ns) >0) {
    # Apply the function
    de <- FUN(Ns, ..., CheckArguments=FALSE)
  } else {
    # or return NA if no points are left
    de <- NA 
  }
  # Return its value according to its class
  if (is.PhyloEntropy(de)) {
    H <- de$Total
  } else {
    H <- de
  }
  return(list(H=H, n=Xr$n))
}
