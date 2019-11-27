####LINDA VALERI####
####09/17/2019######
######3P8122########
####FISHER EXAMPLE##
A <- c(rep(1,8),rep(0,8))
Y <- c(50,62,48,55,58,61,58,56,48,46,54,45,53,46,53,48)

T.obs <- mean(Y[A == 1]) - mean(Y[A == 0])
T.obs

library(ri)
Abold <- genperms(A,maxiter = 12870)

#Abold <- genperms(A)

Abold[, 1:6]


rdist <- rep(NA, times = ncol(Abold))
for (i in 1:ncol(Abold)) {
  A_tilde <- Abold[, i]
  rdist[i] <- mean(Y[A_tilde == 1]) - mean(Y[A_tilde == 0])
}
rdist
hist(rdist)

# p-value
pval <- mean(rdist >= T.obs)
pval
quant <- quantile(rdist,probs = 1-pval)
hist(rdist)
abline(v = quant,col="red")



n.sim = 10000 #number of simulations
T.sim = rep(NA, n.sim)
for (i in 1:n.sim) {
  A.sim = sample(A)
  T.sim[i] = mean(Y[A.sim==1], na.rm=TRUE) -  mean(Y[A.sim==0], na.rm=TRUE)
}
hist(T.sim, xlab = "Test Statistic under Null", main="Randomization Distribution")

#13
mean(T.sim >= T.obs) #approximate p-value
#very close to exact

quant <- quantile(T.sim,probs = 1-pval)
hist(T.sim)
abline(v = quant,col="red")



