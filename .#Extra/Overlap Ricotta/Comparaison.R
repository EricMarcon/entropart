# Comparaison Leps - LC
#======================

# Calculate the overlap, sd=sigma
Overlap <- function(d, sigma) 2*pnorm(d/2, sd=sigma, lower.tail = F)
# Compare the overlap  exp(-ud) where u=1/(sqrt(2)sd)
sigma <- 1/6
d <- seq(from=0, to=max(1, sigma*6), length.out = 100)
plot(Overlap(d, sigma)~d, type="l")
curve(exp(-x/sqrt(2)/sigma), from=0, to=max(1, sigma*6), add=TRUE, col="red")
abline(v=1)
# Minimum similarity
Overlap(1, sigma)


# Comparaison Ricotta - LC
#=========================
# Other presentation: plot LC's Z agaist Leps's Overlap
d0 <- seq(0, 1, .01)
dprime <- 1/3
X <- d0/dprime-(d0/dprime)^2/2+(d0/dprime)^3/6
Y <- 1-exp(-d0/dprime)
plot(Y~X, type ="l")


# Approximation répartition normale
===================================
# OK jusqu'à 2 sd
TailN <- function(x) pnorm(x, lower.tail = F)
Phi1 <- function(x) 1/2 - 1/sqrt(2*pi)*(x)
Phi3 <- function(x) 1/2 - 1/sqrt(2*pi)*(x -x^3/6)
Phi5 <- function(x) 1/2 - 1/sqrt(2*pi)*(x -x^3/6 + x^5/40)
curve(TailN)
curve(Phi1, add=T, col="red")
curve(Phi3, add=T, col="orange")
curve(Phi5, add=T, col="green")

# Grandes distances
TailT1 <- function(x) exp(-x^2)/sqrt(2*pi)*(1/x)
TailT3 <- function(x) exp(-x^2)/sqrt(2*pi)*(1/x -1/x^3)
TailT5 <- function(x) exp(-x^2)/sqrt(2*pi)*(1/x -1/x^3 + 3/x^5)
curve(TailN)
curve(TailT1, add=T, col="red")
curve(TailT3, add=T, col="orange")
curve(TailT5, add=T, col="green")
