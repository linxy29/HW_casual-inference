##############################
#Matching Methods 
##############################

# Install packages
#install.packages("MatchIt")
#install.packages("optmatch")
#install.packages("RItools")
#install.packages("texreg")
#install.packages("Matching")
#install.packages("rgenoud")
#install.packages("Zelig")

# Load packages
library(MatchIt)
library(optmatch)
library(RItools)
library(Matching)
library(rgenoud)
library(texreg)
library(tableone)
library(personalized)
library(Zelig)

## Data
# Load data: a subset of the LaLonde data is loaded with the MatchIt library
data(lalonde)
?lalonde
#To read data into R from outside source, use: read.table, read.csv, etc.
#Or from foreign library, use: read.dta, read.sav, etc.
#Ex: library(foreign)
#Ex: dataname <- read.dta('path/filename.dta')

# Examine data briefly
dim(lalonde)
names(lalonde)
str(lalonde)
head(lalonde)
summary(lalonde)
table(lalonde$treat)

#Check Covariate Balance
vars <- c("age"  ,   "educ" ,   "black"  , "hisp"   , "married" ,"nodegr"  ,"re74"  ,  "re75" )
## Construct a table
covbal0 <- CreateTableOne(vars = vars, strata = "treat", data = lalonde, test = FALSE, smd=TRUE)
print(covbal0, smd = TRUE)

#OTHER FUNCTIONS FOR COVARIATE BALANCE CHECKING
# Check initial balance (from RItools package): 
# How similar are treatment/control groups on X?
#?xBalance
#xBalance(treat ~ . - (re78 + treat), data =lalonde, report=c("adj.means","std.diffs")) 
#xBalance(treat ~ . - (re78 + treat), data =lalonde, report=c("all")) 
#?balanceUV
#RE75.BAL  <- balanceUV(lalonde$re75[lalonde$treat==1],lalonde$re75[lalonde$treat!=1])
#summary(RE75.BAL)



## Example USING MATCHIT PACKAGE 

# 1. Implement nearest neighbor matching on the propensity score (NN match on PS score: 1:1, greedy, without replacement or calipers)
psmatch1 <- matchit(treat ~ age + educ + black + hisp + nodegr + married + re74 + re75, 
                    distance="logit", method = "nearest", discard = "control", data = lalonde)

# 2. Check balance
summary(psmatch1, standardize=TRUE)
#compare with initial sample
print(covbal0, smd = TRUE)

plot(psmatch1)
plot(psmatch1, type="hist")
par(mfrow = c(1, 1))
plot(psmatch1, type="jitter", interactive=FALSE)


## More on nearest neighbor matching and the propensity score

# Estimating a propensity score, implementing a matching procedure, checking balance

# Optimal (instead of greedy) matching on PS score: 1:1, without replacement or calipers
# Now Add polynomials and interactions to propensity score model
psmatch2 <- matchit(treat ~ age + educ + black + hisp + nodegr + married + re74 +re75
                    + I(age^2) + nodegr:re74,  distance = "logit", method = "optimal", data = lalonde)
#authors of the R package say to ignore this warning..

# Check balance
summary(psmatch2, standardize=TRUE)
plot(psmatch2)
plot(psmatch2, type="jitter", interactive = FALSE)
plot(psmatch2, type="hist")

# 3. Create matched data for analysis
psmatch1.data <- match.data(psmatch1)
psmatch2.data <- match.data(psmatch2)


# 4. Analyize: Difference of means
# Outcome Regression Analysis: AVERAGE TREATMENT EFFECT ON THE TREATED
psmatch1.mod <- lm(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75, 
                   data = psmatch1.data)
summary(psmatch1.mod)  #beta of treat：ATT, among the treated
?zelig
psmatch1.mod.zelig <- zelig(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75,  data = psmatch1.data,  model = "ls")
summary(psmatch1.mod.zelig)

psmatch1.mod <- lm(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75 + 
                     I(age^2) + I(educ^2) + re74 + nodegr:re74, 
                   data = psmatch1.data)
summary(psmatch1.mod) 
psmatch2.mod <- lm(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75 + 
                     I(age^2) + I(educ^2) + re74 + nodegr:re74, 
                   data = psmatch2.data) # optimal matching
summary(psmatch2.mod)


### AVERAGE TREATMENT EFFECT 
z.out1 <- zelig(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75 + 
                 I(age^2) + I(educ^2) + re74 + nodegr:re74, data = match.data(psmatch1), model = "ls")
x.out1 <- setx(z.out1, treat=0)
x1.out1 <- setx(z.out1, treat=1)
s.out1 <- Zelig:::sim(z.out1, x = x.out1, x1 = x1.out1)
summary(s.out1)

z.out2 <- zelig(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75 + 
                  I(age^2) + I(educ^2) + re74 + nodegr:re74, data = match.data(psmatch2), model = "ls")
x.out2 <- setx(z.out2, treat=0)
x1.out2 <- setx(z.out2, treat=1)
s.out2 <- Zelig:::sim(z.out2, x = x.out2, x1 = x1.out2)
summary(s.out2)

#COMPARE WITH ASSOCIATION
# Treatment effect without matching: difference of means
association <- lm(re78 ~ treat, data=lalonde)
summary(association)

#COMPARE WITH DIRECT REGRESSION ADJUSTMENT
# Treatment effect without matching: regression
unmatched <- lm(re78 ~ treat + age +I(age^2)+ educ + I(educ^2) + black + hisp + nodegr + married +nodegr*re74+ re75, 
                data = lalonde)
summary(unmatched)


## Final comparison
association$coeff["treat"]
unmatched$coeff["treat"]
summary(s.out1)
summary(s.out2)



####END OF LAB####









####OPTIMAL MATCHING WITH GENMATCH


#This function finds optimal balance using multivariate matching where a search 
#algorithm determines the weight each covariate is given. Balance is determined by 
#examining cumulative probability distribution functions of a variety of standardized 
#statistics. By default, these statistics include t-tests and Kolmogorov-Smirnov tests.
#A variety of descriptive statistics based on empirical- QQ (eQQ) plots can also be used 
#or any user provided measure of balance. The statistics are not used to conduct formal
#hypothesis tests, because no measure of balance is a monotonic function of bias and 
#because balance should be maximized without limit. 
#The object returned by GenMatch can be supplied to the Match function 
#(via the Weight.matrix option) to obtain causal estimates.

?GenMatch
attach(lalonde)
#The covariates we want to match on
X = cbind(age, educ, black, hisp, married, nodegr, re75, re74)
#The covariates we want to obtain balance on 
BalanceMat <- cbind(age,I(age^2), educ,I(educ^2), black, hisp, married, nodegr, re75, re74, I(nodegr*re74))
#
#Let's call GenMatch() to find the optimal weight to give each
#covariate in 'X' so as we have achieved balance on the covariates in
#'BalanceMat'. This is only an example so we want GenMatch to be quick
#so the population size has been set to be only 16 via the 'pop.size'
#option. This is *WAY* too small for actual problems.
#For details see http://sekhon.berkeley.edu/papers/MatchingJSS.pdf.
#
out <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat, estimand="ATE", M=1,
                   pop.size=200, max.generations=10, wait.generations=1)
#The outcome variable
Y=re78/1000
#


#
#Let's determine if balance has actually been obtained on the variables of interest
#
mb <- MatchBalance(treat~age +educ+black+ hisp+ married+ nodegr+ re75+ re74+ I(re74*re75) + I(re74*nodegr),
                   match.out=out, nboots=500)


# Now that GenMatch() has found the optimal weights, let's estimate
# our causal effect of interest using those weights
#
mout <- Match(Y=Y, Tr=treat, X=X, estimand="ATE", Weight.matrix=out)

psmatch3.mod <-mout

summary(psmatch3.mod)


#########OTHER MATCHING TECHNIQUES WITH MATCHIT

# NN match on PS score: k:1, greedy, without replacement or calipers
psmatch4 <- matchit(treat ~ age + educ + black + hisp + nodegr + married + re75
                    + I(age^2) + nodegr:re75,  distance = "logit", method = "nearest", discard= "control", 
                    ratio = 3, data = lalonde)

# Check balance
summary(psmatch4, standardize=TRUE)
plot(psmatch4)
plot(psmatch4, type="hist")
par(mfrow = c(1, 1))
plot(psmatch4, type="jitter", interactive = FALSE)

# Create matched data for analysis
psmatch4.data <- match.data(psmatch4)

# NN match on PS score: k:1, with replacement (greedy doesn't matter), no calipers
psmatch5 <- matchit(treat ~ age + educ + black + hisp + nodegr + married + re75+ I(age^2) + nodegr:re75, 
                    distance = "logit", method = "nearest", discard= "control", 
                    ratio = 3, replace= TRUE, data = lalonde)
# Check balance
summary(psmatch5, standardize=TRUE)
plot(psmatch5)
plot(psmatch5, type="hist")
par(mfrow = c(1, 1))
plot(psmatch5, type="jitter", interactive = FALSE)

#now ratio 1
psmatch6 <- matchit(treat ~ age + I(age^2) + educ + I(educ^2) + black +
                      hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
                      u74 + u75,     ratio = 1, replace= TRUE, data = lalonde)

# Check balance
summary(psmatch6, standardize=TRUE)
plot(psmatch6)
plot(psmatch6, type="hist")
par(mfrow = c(1, 1))
plot(psmatch6, type="jitter", interactive = FALSE)

# NN match on PS score: k:1, with replacement, add calipers (instead of discarding outside common support)
psmatch7 <- matchit(treat ~ age + educ + black + hisp + nodegr + married + re75 + 
                      I(age^2) + nodegr:re75, 
                    distance = "logit", method = "nearest", caliper = 0.3, 
                    ratio = 3, replace= TRUE, data = lalonde)

# Check balance
summary(psmatch7, standardize=TRUE)
plot(psmatch7, type="jitter", interactive = FALSE)
plot(psmatch7, type="hist")

# Create matched data for analysis
psmatch4.data <- match.data(psmatch4)
psmatch5.data <- match.data(psmatch5)
psmatch6.data <- match.data(psmatch6)
psmatch7.data <- match.data(psmatch7)

#OUTCOME ANALYSIS
psmatch4.mod <- lm(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75 + 
                     I(age^2) + I(educ^2) + educ:re74, 
                   data = psmatch4.data)
summary(psmatch4.mod)

psmatch5.mod <- lm(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75 + 
                     I(age^2) + I(educ^2) + educ:re74, 
                   data = psmatch5.data, weights = weights)
summary(psmatch5.mod)
# Added weights option, weights = weights
# matchit automatically created frequency weights (named `weights') in the matched data frame
psmatch6.mod <- lm(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75 + 
                     I(age^2) + I(educ^2) + educ:re74, 
                   data = psmatch5.data, weights = weights)
summary(psmatch6.mod)


psmatch6.mod <- lm(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75 + 
                     I(age^2) + I(educ^2) + educ:re74, 
                   data = psmatch6.data, weights = weights)
summary(psmatch6.mod)

psmatch7.mod <- lm(re78 ~ treat + age + educ + black + hisp + nodegr + married + re74 + re75 + 
                     I(age^2) + I(educ^2) + educ:re74, 
                   data = psmatch7.data, weights = weights)
summary(psmatch7.mod)



