rCommunity <- 
function(n, size = sum(NorP), NorP = 1, BootstrapMethod = "Chao2015",
         S = 300, Distribution = "lnorm", sd = 1, prob = 0.1, alpha=40, 
         CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  Ps <- Ns <- NULL
  # Draw probabilities, except for logseries: draw abundances.
  if (length(NorP) == 1) {
    # Draw in a distribution
    if (Distribution == "lseries") {
      # Abundance of a species in a logseries distribution of size N and parameter alpha
      LSabundance <- function(N, alpha) {
        # adapted from Dan Lunn, http://www.stats.ox.ac.uk/~dlunn/BS1_05/BS1_Rcode.pdf
        # Fisher's x is log-series 1-theta
        x <- N/(N+alpha)
        # Draw a random number between 0 and 1
        u <- stats::runif(1)
        # k is the number of individuals to draw
        k <- 1
        # Calculate the probability at k=1
        P <- -x/log(1-x)
        # Store it in the distribution function
        F <- P
        # Repeat while the cumulated probabiilty is below u
        while (u >= F) {
          # Probability at k+1 obtained from that at k
          P <- P*k*x/(k+1)
          # Increment k
          k <- k+1
          # Increment the cumulated probability
          F <- F+P
        }
        return(k)
      }
      # Draw the abundances of the number of species corresponding to size and alpha
      Ns <-replicate(n, replicate(round(-alpha*log(alpha/(size+alpha))), LSabundance(size, alpha)))
    } else {
      # Other distributions: draw probabilities
      Ps <- switch(Distribution,
                   geom = prob/(1-(1-prob)^S)*(1-prob)^(0:(S-1)),
                   lnorm = (stats::rlnorm(S, 0, sd) -> Nslnorm)/sum(Nslnorm),
                   bstick = c(cuts <- sort(stats::runif(S-1)), 1)- c(0, cuts)
                  )
    }
    
  } else {
    # Subsample
    if (abs(sum(NorP) - 1) < length(NorP)*.Machine$double.eps) {
      # Probabilities sum to 1, allowing rounding error
      Ps <- NorP    
    } else {
      # Abundances: Generate Ps according to the chosen method
      if (BootstrapMethod == "Chao2015") {
        Ps <- as.ProbaVector(NorP, Correction = "Chao2015", Unveiling = "geom", CheckArguments = FALSE)
      }
      if (BootstrapMethod == "Chao2013") {
        Ps <- as.ProbaVector(NorP, Correction = "Chao2013", Unveiling = "unif", CheckArguments = FALSE)
      }
      if (BootstrapMethod == "Marcon") {
        Ps <- NorP/sum(NorP)
      }
    }
  }
  
  if (is.null(Ps) & is.null(Ns)) {
    warning ("The distribution to simulate has not been recognized")
    return(NA)
  }
  
  # Generate communities according to Ps
  if (is.null(Ns)) {
    # Draw a multinomial sample from Ps except if Ns has already been obtained (e.g.: lseries)
    Ns <- stats::rmultinom(n, size, Ps)
  }
 
  if (n > 1) {
    # Return a MetaCommunity
    return(MetaCommunity(Ns))      
  } else {
    # Return a vector if a single community has been simulated
    return(as.AbdVector(Ns, Round = TRUE))
  }
}
