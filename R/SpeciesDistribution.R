as.SpeciesDistribution <-
function (x, ...)
{
  UseMethod("as.SpeciesDistribution")
}


as.SpeciesDistribution.data.frame <-
function (x, ...) 
{
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")

  # Try to save the names before applying as.numeric
  spNames <- names(x)
  spD <- as.numeric(as.matrix(x))
  if (length(spNames) == length(spD))
    names(spD) <- spNames
  
  class(spD) <- c("SpeciesDistribution", class(spD))
  return(spD)
}


as.SpeciesDistribution.numeric <-
function (x, ...) 
{
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")
  
  spD <- x
  class(spD) <- c("SpeciesDistribution", class(spD))
  return(spD)
}


as.SpeciesDistribution.integer <-
function (x, ...) 
{
  return(as.SpeciesDistribution.numeric(x))
}


is.SpeciesDistribution <-
function (x)
{
  inherits(x, "SpeciesDistribution")
}


as.ProbaVector <-
function (x, ...) 
{
  UseMethod("as.ProbaVector")
}


as.ProbaVector.data.frame  <-
function (x, ...) 
{
  spD <- as.SpeciesDistribution(x, ...)
  
  return(as.ProbaVector.numeric(spD, CheckArguments=FALSE))
}


# Utilities for as.ProbaVector.numeric. Not exported. ####
# Solve the theta parameter of Chao et al. (2015)
theta_solve <- function(theta, Ps, Ns, N, C, CD2){
  # Code inspired from JADE function DetAbu(), http://esapubs.org/archive/ecol/E096/107/JADE.R
  lambda <- (1-C) / sum(Ps * exp(-theta*Ns))
  return(abs(sum((Ps * (1 - lambda * exp(-theta*Ns)))^2) - sum(choose(Ns,2)/choose(N,2)) + CD2))
}
# Solve the beta parameter of Chao et al. (2015)
beta_solve <- function(beta, r, i){
  # Code inspired from JADE function UndAbu(), http://esapubs.org/archive/ecol/E096/107/JADE.R
  return(abs(sum(beta^i)^2 / sum((beta^i)^2) - r))
}
# Unobserved species distribution
estimate_Ps0 <- function(Unveiling, PsTuned, S0, C, CD2){
  Ps0 <- NA
  if (Unveiling == "geom") {
    if (S0 == 1) {
      # A single unobserved species
      Ps0 <- 1-C
    } else {
      r <- (1-C)^2/CD2
      i <- 1:S0
      beta <-  tryCatch(stats::optimize(beta_solve, lower=(r-1)/(r+1), upper=1, tol=.Machine$double.eps, r, i)$min, 
                        error = function(e) {(r-1)/(r+1)})
      alpha <- (1-C) / sum(beta^i)
      Ps0 <- alpha * beta^i
      # Sometimes fails when the distribution is very uneven (sometimes r < 1) 
      # Then, go back to the uniform distribution
      if (any(is.na(Ps0)) | any(Ps0 <= 0)) Unveiling <- "unif"
    }
  }      
  if (Unveiling == "unif") {
    # Add S0 unobserved species with equal probabilities
    Ps0 <- rep((1-sum(PsTuned))/S0, S0)
  }
  if (any(is.na(Ps0))) {
    warning("Unveiling method was not recognized")
    return(NA)
  } else {
    names(Ps0) <- paste("UnobsSp", 1:(length(Ps0)), sep="")
    return(Ps0)
  }         
}
# Rarefaction bias
rarefaction_bias <- function(S0, Ns, PsTuned, C, CD2, q, Unveiling, Target) {
  Ns <- Ns[Ns>0]
  N <- sum(Ns)
  # Unobserved species
  Ps0 <- estimate_Ps0(Unveiling, PsTuned, S0, C, CD2)
  # Full distribution of probabilities
  Ps <- c(PsTuned, Ps0)
  # AbdFreqCount at Level = N
  Sn <- sapply(1:N, function(nu) sum(exp(lchoose(N, nu) + nu*log(Ps) + (N-nu)*log(1-Ps))))
  # Get Entropy at Level=N and calculate the bias
  if (q == 1) {
    Bias <- abs(sum(-(1:N)/N * log((1:N)/N) * Sn) - Target)
  } else {
    Bias <- abs((sum(((1:N)/N)^q * Sn) - 1) / (1-q) - Target)
  }
  return(Bias)
}
# end of utilities ####

as.ProbaVector.numeric <-
function (x, Correction = "None", Unveiling = "None", RCorrection = "Jackknife", 
          JackOver = FALSE, CEstimator = "ZhangHuang", q = 0, ..., CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Try to save the names before applying as.vector
  spNames <- names(x)
  spD <- as.SpeciesDistribution(as.vector(x), ...)
  if (length(spNames) == length(spD))
    names(spD) <- spNames

  if (Correction == "None") { 
    spD <- spD/sum(spD)
  } else {
    # Integer abundances are required
    if (!is.IntValues(spD)) warning("Integer abundance values are required to estimate community probabilities. Abundances have been rounded.")
    NsInt <- round(spD)

    # Eliminate 0 and calculate elementary statistics
    Ns <- NsInt[NsInt > 0]
    S <- length(Ns)
    N <- sum(Ns)
    Ps <- Ns/N
    # Sample coverage
    C <- Coverage(Ns, Estimator=CEstimator)
    if (Correction == "Chao2015" | Unveiling == "Chao2015" | RCorrection == "Rarefy") {
      # Sample coverage of order 2 required
      Singletons <- sum(Ns==1)
      Doubletons <- sum(Ns==2)
      if (Doubletons==0) {
        Singletons <- max(Singletons - 1, 0)
        Doubletons <- 1
      }
      Tripletons <- max(sum(Ns==3), 1)
      # 1 minus sample coverage (i.e. Coverage Deficit) of order 2
      CD2 <- Doubletons / choose(N, 2) * ((N-2)*Doubletons / ((N-2)*Doubletons + 3*Tripletons))^2
    }
    
    # Tune the probabilities of observed species
    if (C == 0 | C == 1) {
      # Sample coverage equal to 1: do not tune. If 0, unable to tune.
      PsTuned <- Ps
    } else {
      PsTuned <- NA
      if (Correction == "ChaoShen") {
        PsTuned <- C*Ps
      }
      if (Correction == "Chao2013") {
        # Single parameter estimation, Chao et al. (2013)
        denominator <- sum(Ps*(1-Ps)^N)
        if (denominator == 0) {
          # N is too big so denominator equals 0. Just multiply by C.
          PsTuned <- Ps*C
        } else {
          # General case
          lambda <- (1 - C)/denominator
          PsTuned <- Ps*(1 - lambda*(1-Ps)^N)
        }      
      } 
      if (Correction == "Chao2015")  {
        # Two parameters, Chao et al. (2015). 
        # Estimate theta. Set it to 1 if impossible
        theta <- tryCatch(stats::optimize(theta_solve, interval=c(0,1), Ps, Ns, N, C, CD2)$min, 
                          error = function(e) {1})
        lambda <- (1-C) / sum(Ps * exp(-theta*Ns))
        PsTuned <- Ps * (1 - lambda * exp(-theta*Ns))
      }
      if (any(is.na(PsTuned))) {
        warning("Correction was not recognized")
        return (NA)
      }
    }
    names(PsTuned) <- names(spD[spD > 0])
    
    # Estimate the number of unobserved species
    if (RCorrection == "Rarefy") {
      if (Unveiling == "None")
        stop("Arguments RCorrection='Rarefy' and Unveiling='None' are not compatible")
      # Estimation of the number of unobserved species to initialize optimization
      S0 <- bcRichness(Ns, Correction="Jackknife") - S
      # Estimate the number of unobserved species by iterations
      Target <- Tsallis(Ns, q=q, Correction="None", CheckArguments = FALSE)
      S0 <- round(tryCatch(stats::optimize(rarefaction_bias, interval=c(0, 2*S0), Ns, PsTuned, C, CD2, q, Unveiling, Target)$minimum,
                           error = function(e) {S0}))
    } else {
      Sestimate <- ceiling(bcRichness(Ns, Correction=RCorrection, JackOver=JackOver))
      S0 <- Sestimate - S
    }
    
    # Distribution of unobserved species
    if (S0) {
      if (Unveiling == "None") {
        spD <- PsTuned
      } else {
        spD <- c(PsTuned, estimate_Ps0(Unveiling, PsTuned, S0, C, CD2))
      }
    } else {
      spD <- PsTuned
    }
    spD <- as.SpeciesDistribution(spD, ...)
  }
  class(spD) <- c("ProbaVector", class(spD))
  return(spD)
}


as.ProbaVector.integer <-
function (x, Correction = "None", Unveiling = "None", RCorrection = "Jackknife", 
          JackOver = FALSE, CEstimator = "ZhangHuang", q = 0, ..., CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return(as.ProbaVector.numeric(x, Correction=Correction, Unveiling=Unveiling, RCorrection=RCorrection, 
                                JackOver=JackOver, CEstimator=CEstimator, q=q, ..., CheckArguments=FALSE))
}


is.ProbaVector <-
function (x) 
{
  inherits(x, "ProbaVector")
}


as.AbdVector <-
function (x, ...)
{
  UseMethod("as.AbdVector")
}


as.AbdVector.data.frame <-
function (x, Round = TRUE, ...) 
{
  # Try to save the names before applying as.vector
  spNames <- names(x)
  
  if (Round) {
    intx <- as.integer(as.matrix(round(x)))
    spD <- as.SpeciesDistribution(as.vector(intx), ...)
  } else {      
    spD <- as.SpeciesDistribution(as.vector(x), ...)
  }

  # Restore the names
  if (length(spNames) == length(spD))
    names(spD) <- spNames
  
  class(spD) <- c("AbdVector", class(spD))
  return(spD)
}


as.AbdVector.numeric <-
function (x, Round = TRUE, ...) 
{
  # Try to save the names before applyinf as.vector
  spNames <- names(x)

  if (Round) {
    intx <- as.integer(round(x))
    spD <- as.SpeciesDistribution(as.vector(intx), ...)
  } else {      
    spD <- as.SpeciesDistribution(as.vector(x), ...)
  }
  
  # Restore the names
  if (length(spNames) == length(spD))
    names(spD) <- spNames
  
  class(spD) <- c("AbdVector", class(spD))
  return(spD)
}


as.AbdVector.integer <-
function (x, ...) 
{
  return(as.AbdVector.numeric(x))
}


is.AbdVector <-
function (x) 
{
  inherits(x, "AbdVector")
}


plot.SpeciesDistribution <-
function(x, ..., Distribution = NULL, 
         type = "b", log = "y", main = NULL, xlab = "Rank", ylab = NULL) 
{
  # Eliminate zeros and sort
  Ns <- sort(x[x > 0], decreasing = TRUE)
  N <- sum(Ns)
  S <- length(Ns)
  
  # Prepare ylab
  if (is.null(ylab)) {
    if (is.ProbaVector(x)) {
      ylab <- "Probability"
    } else {
      ylab <- "Abundance" 
    }
  }
  
  graphics::plot(Ns, type=type, log=log, main=main, xlab=xlab, ylab=ylab, axes=FALSE, ...)
  # x axis ticks must start from 1
  graphics::axis(1, graphics::axTicks(1)+1)
  graphics::axis(2)
  graphics::box()
  
  if (!is.null(Distribution)) {
    if (Distribution == "lnorm") {
      FittedRAC <- RAClnorm(Ns)
      graphics::lines(FittedRAC$Rank, FittedRAC$Abundance, col = "red")
      return(list(mu = FittedRAC$mu, sigma = FittedRAC$sigma))
    }
    if (Distribution == "geom") {
      FittedRAC <- RACgeom(Ns)
      graphics::lines(FittedRAC$Rank, FittedRAC$Abundance, col = "red")
      return(list(prob = FittedRAC$prob))
    }
    if (Distribution == "lseries") {
      FittedRAC <- RACgeom(Ns)
      graphics::lines(FittedRAC$Rank, FittedRAC$Abundance, col = "red")
      return(list(alpha = FittedRAC$alpha))
    }
    if (Distribution == "bstick") {
      FittedRAC <- RACgeom(Ns)
      graphics::lines(FittedRAC$Rank, FittedRAC$Abundance, col = "red")
      return(list(max = FittedRAC$max))
    }
    warning("The distribution to fit has not been recognized")
    return(NA)
  }
}


autoplot.SpeciesDistribution <-
function(object, ..., Distribution = NULL, 
         ylog = TRUE, main = NULL, xlab = "Rank", ylab = NULL) 
{
  # Eliminate zeros and sort
  Ns <- sort(object[object > 0], decreasing = TRUE)
  N <- sum(Ns)
  S <- length(Ns)
  
  # Transform data into df
  df <- data.frame(Rank=1:S, Ns)
  
  # Prepare ylab
  if (is.null(ylab)) {
    if (is.ProbaVector(object)) {
      ylab <- "Probability"
    } else {
      ylab <- "Abundance" 
    }
  }

  # Plot. X-axis starts at 0.01 to avoid the 0 X-label.
  thePlot <- ggplot2::ggplot() +
    ggplot2::geom_point(data=df, mapping=ggplot2::aes_(x=~Rank, y=~Ns)) +
    ggplot2::scale_x_continuous(limits=c(0.01, S), expand=c(0, 0)) +
    ggplot2::labs(title=main, x=xlab, y=ylab)
  
  # Log Y-axis
  if (ylog) thePlot <- thePlot + ggplot2::scale_y_log10() 
  
  OK <- TRUE
  # Fit distributions
  if (!is.null(Distribution)) {
    OK <- FALSE
    if (Distribution == "lnorm") {
      FittedRAC <- RAClnorm(Ns)
      OK <- TRUE
    }
    if (Distribution == "geom") {
      FittedRAC <- RACgeom(Ns)
      OK <- TRUE
    }
    if (Distribution == "lseries") {
      FittedRAC <- RAClseries(Ns)
      OK <- TRUE
    }
    if (Distribution == "bstick") {
      FittedRAC <- RACbstick(Ns)
      OK <- TRUE
    }
    if (OK) {
      # Add the adjusted curve to the plot
      thePlot <- thePlot + 
        ggplot2::geom_line(data=with(FittedRAC, data.frame(Rank, Abundance)), mapping=ggplot2::aes_(x=~Rank, y=~Abundance, col="red")) + 
        ggplot2::theme(legend.position="none")
      # Add fitted parameters to the attributes of the plot
      if (Distribution == "lnorm") {
        attr(thePlot, "mu") <- FittedRAC$mu
        attr(thePlot, "sigma") <- FittedRAC$sigma
      }
      if (Distribution == "geom") {
        attr(thePlot, "prob") <- FittedRAC$prob
      }
      if (Distribution == "lseries") {
        attr(thePlot, "alpha") <- FittedRAC$alpha
      }
      if (Distribution == "bstick") {
        attr(thePlot, "max") <- FittedRAC$max
      }
    }
  } 
  if (OK) {
    # Return the plot
    return(thePlot)
  } else {
    # Return NA with a warning
    warning("The distribution to fit has not been recognized")
    return(NA)   
  }
}


is.IntValues <-
function (Ns)
{
  NsInt <- round(Ns)
  # Return TRUE if no value in Ns has been modified by rounding
  return(!(any(abs(NsInt-Ns) > sum(Ns)*.Machine$double.eps)))
}


RAClnorm <- 
function (Ns, CheckArguments = TRUE)
{
  # Eliminate zeros
  Ns <- sort(Ns[Ns > 0], decreasing = TRUE)
  S <- length(Ns)
  
  # Fit a lognormal distribution
  mu <- mean(log(Ns))
  sigma <- stats::sd(log(Ns))
  # Unique values
  Ns1 <- unique(Ns)
  Rank <- S*(1-stats::pnorm(log(Ns1), mu, sigma))
  
  return(list(Rank=Rank, Abundance=Ns1, mu=mu, sigma=sigma))
}


RACgeom <- 
function (Ns, CheckArguments = TRUE)
{
  # Eliminate zeros
  Ns <- sort(Ns[Ns > 0], decreasing = TRUE)
  S <- length(Ns)
  
  # Fit a geometric distribution
  lNs <- log(Ns)
  Rank <- 1:S
  reg <- stats::lm(lNs~Rank)
  
  return(list(Rank=Rank, Abundance=exp(reg$coefficients[1]+reg$coefficients[2]*Rank), prob=as.numeric(-reg$coefficients[2])))
}


RAClseries <- 
function (Ns, CheckArguments = TRUE)
{
  # Eliminate zeros
  Ns <- sort(Ns[Ns > 0], decreasing = TRUE)
  N <- sum(Ns)
  
  # Evaluate alpha
  alpha <- vegan::fisher.alpha(Ns)
  # May (1975) Ecology and evolution of communities, Harvard University Press
  sei <- function(t) exp(-t)/t
  Rank <- vapply(unique(Ns), function(x) {
    n <- x * log(1 + alpha/N)
    f <- stats::integrate(sei, n, Inf)
    fv <- f[["value"]]
    return(alpha * fv)}
    , 0)
  
  return(list(Rank=Rank, Abundance=unique(Ns), alpha=alpha))
}


RACbstick <- 
function (Ns, CheckArguments = TRUE)
{
  # Eliminate zeros
  Ns <- sort(Ns[Ns > 0], decreasing = TRUE)
  N <- sum(Ns)
  S <- length(Ns)
  
  # Fit a broken stick
  f1 <- sort(cumsum(1/(S:1)), decreasing = TRUE)
  Rank <- 1:S
  Abundance <- N*f1/sum(f1)
  
  return(list(Rank=Rank, Abundance=Abundance, max=max(Abundance)))
}
