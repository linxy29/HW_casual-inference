#####MSM EXAMPLE CODE FOR A SINGLE EXPOSURE TIME POINT#######
###################USING GARDASIL DATA#######################
######################LINDA VALERI###########################

####NOTE: FOR HOMEWORK NEED TO EXPAND THIS FOR MULTIPLE EXPOSURE TIME POINTS####

x <- read.table("./gardasil.dat.txt",header = T)
head(x)

#Y: Completed
#A: PracticeType

x$PracticeType <- as.factor(x$PracticeType)
summary(x$PracticeType)
#0 pediatrics
#1 family practice 
#2 OB-GYN



#C: Age, Race, InsuranceType, MedAssist, LocationType

#Propensity score
x$PracticeType_bin <- as.numeric(as.numeric(x$PracticeType)>1)
ps.model<-glm(PracticeType_bin~Age +  as.factor(Race) +   as.factor(InsuranceType) + as.factor(Location) ,data=x, family = binomial)
summary(ps.model)


# install R package for visualization of overlap

library(personalized)
prop.func <- function(x, trt)
{
  # fit propensity score model
  propens.model <- glm(trt~Age +  as.factor(Race) +   as.factor(InsuranceType) + as.factor(Location) ,data=x, family = binomial)
  pi.x <- predict(propens.model, type = "response")
  pi.x
}

check.overlap(x = x,
              trt = x$PracticeType_bin,
              propensity.func = prop.func)


# now add density plot with histogram
check.overlap(x = x,
              trt = x$PracticeType_bin,
              type = "both",
              propensity.func = prop.func)




library(tableone)
x$InsuranceType <-as.factor(x$InsuranceType)
x$Location <-as.factor(x$Location)
x$Race <-as.factor(x$Race)
vars <- c("Age" , "Race", "InsuranceType" ,"Location")

## Construct a table
tab <- CreateTableOne(vars = vars, strata = "PracticeType_bin", data = x, test = FALSE)


## Show table with SMD

print(tab, smd = TRUE)

summary(x$Age)
dim(x)
dim(x[which(x$Age <=18),])
dim(x[which(x$PracticeType !=0),])
summary(x$Age[which(x$PracticeType !=0)])
summary(x$Age[which(x$PracticeType ==0)])

#################REGRESSION ADJUSTMENT###################

reg_adj <- glm(Completed~PracticeType_bin+ Age +as.factor(Race) +as.factor(InsuranceType) + as.factor(Location),data=x, family = binomial)
summary(reg_adj)

#################MARGINAL STRUCTURAL MODEL##############
nboot <- 1000
n <- nrow(x)
library(survey)
b.holder <- rep(NA)
for (b in 1:nboot){
  set.seed(123+b)
  S.b <- sample(1:n, size = n, replace = TRUE)
  data.b <- x[S.b, ]
  ps.model.b <- glm(PracticeType_bin~Age +  as.factor(Race) +   as.factor(InsuranceType) + as.factor(Location) ,data=data.b, family = binomial)

ps <- predict(ps.model.b, type="response") #gets the propensity scores for each unit, based on the model
data.b$est.w <- ifelse(data.b$PracticeType_bin == 1, 1/ps, 1/(1 - ps))
msm <-svyglm(Completed~PracticeType_bin ,family=binomial(link = "logit"),design = svydesign(~ 1, weights = ~ est.w,data=data.b))
b.holder[b] <- msm$coef[2]
}

quantile(b.holder, probs = c(0.5,0.025,0.975))


#############ASSOCIATION#############

reg <- glm(Completed~PracticeType_bin,data=x, family = binomial)
summary(reg)


#point estimates are different. Final conclusion similar.

