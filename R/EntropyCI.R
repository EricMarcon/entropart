EntropyCI <-
function(FUN, Simulations = 100, Ns, BootstrapMethod = "Chao2015", ShowProgressBar = TRUE, ..., CheckArguments = TRUE) 
{
  if (CheckArguments) {
    CheckentropartArguments()
  }
  
  RealEst <- FUN(Ns, ..., CheckArguments = FALSE)
  
  # SimulateEntropy resamples the community and calculates FUN
  SimulateEntropy <- function(Progress) {
    # Draw Ns from a multinomial distribution
    SimNs <- rCommunity(1, size=sum(Ns), NorP=Ns, BootstrapMethod=BootstrapMethod, CheckArguments = FALSE)
    # FUN(simulated data)
    NewSim <- FUN(SimNs, ..., CheckArguments = FALSE)
    # update progress bar
    if(ShowProgressBar & interactive()) 
      utils::setTxtProgressBar(ProgressBar, Progress)
    return(NewSim)
  }
  
  # Simulate entropy
  ProgressBar <- utils::txtProgressBar(min=0, max=Simulations)
  # Simulated values
  RawSimulatedEntropy <- sapply(1:Simulations, SimulateEntropy)
  close(ProgressBar)
  
  # Recenter entropy
  return(RawSimulatedEntropy + RealEst - mean(RawSimulatedEntropy))
}
