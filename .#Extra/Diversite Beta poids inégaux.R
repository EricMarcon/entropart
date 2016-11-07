effectifs <- c(5,2,3,2,4, 0,0,0,0,    0,0,0,0,0, 7,2 ,0,0,    0,0,0,0,0,0,0, 2,2)
abondance <- matrix(effectifs, ncol=3)
poids <- c(1.1,1,1)
q.seq = seq(0, 2, 0.1)

MC <- MetaCommunity(abondance, poids)
cbind(MC$Psi, MC$Ps)
MC$Wi
plot(Profile <- DivProfile(q.seq, MC, Biased = TRUE))
plot(CommunityProfile(Diversity, poids/sum(poids), q.seq),  type="l")
plot(y=Profile$TotalBetaDiversity, 
     x=CommunityProfile(Diversity, poids/sum(poids), q.seq)$y, 
     type="l")
abline(a=0, b=1, col="red")
DivPart(0, MC, Biased=TRUE)$CommunityAlphaDiversities