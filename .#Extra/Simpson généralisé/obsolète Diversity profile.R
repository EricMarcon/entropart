# Average number of singletons from a uniform distribution (sample size, number of species)
NumSing <- function(v, K) {
  Ns1000 <- rmultinom(n=1000, size=v, prob=rep(1/K, K))
  NumSing <-mean(apply(Ns1000, 2, Singletons <- function(Ns) sum(Ns==1)))
}

# Sample size from 1 to 1000, 100 or 25 species
v.seq <- 1:1000
NS100 <- sapply(v.seq, NumSing, 100)
plot(NS100 ~ v.seq, type="l", xlab="v", ylab="Singletons")
NS25 <- sapply(v.seq, NumSing, 25)
lines(NS25 ~ v.seq, lty=2)

# Sample coverage
Coverage100 <- 1-atan(NS100/v.seq)
plot(Coverage100 ~ v.seq, type="l", xlab="v", ylab="Coverage")
Coverage25 <- 1-atan(NS25/v.seq)
lines(Coverage25 ~ v.seq, lty=2)

# Numbers equivalent
plot(1/(1-(NS100/(v.seq+1))^(1/v.seq)) ~ v.seq, type="l")
