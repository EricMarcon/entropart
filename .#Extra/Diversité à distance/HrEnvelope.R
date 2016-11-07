# Envelope of simulated entropy along a vector of distances
HrEnvelope <-  function(X, r, nPoints, FUN, ..., NumberOfSimulations = 100, Alpha = 0.05, 
           SimulationType = "RandomLabeling", Global = FALSE, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()  

  # Choose the null hypothesis
  SimulatedPP <- switch (SimulationType,
                         RandomPosition = expression(rRandomPositionK(X, CheckArguments = FALSE)),
                         RandomLabeling = expression(rRandomLabeling(X, CheckArguments = FALSE))
  )
  if (is.null(SimulatedPP))
    stop(paste("The null hypothesis", sQuote(SimulationType), "has not been recognized."))

  # local envelope, keep extreme values for lo and hi (nrank=1)
  Envelope <- envelope(X, fun=Hr, nsim=NumberOfSimulations, nrank=1,
                       r=r, nPoints=nPoints, FUN=FUN, ...,
                       CheckArguments = FALSE,
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
