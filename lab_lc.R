######MEDIATION ANALYSIS LAB#######

x<-read.csv("./sim_data_lung.csv")
p<-mean(x$case)
pi<-0.082642
x$w <- ifelse(x$case==1,pi/p,(1-pi)/(1-p)) #weight for y=1 and y=0#

library(mediation)
library(survey)
yreg<-svyglm(case ~ snp*smoking+  sex + colgrad +age, family=binomial(link="logit"),design = svydesign(~ 1, weights = ~ w,data=x))
summary(yreg)
mreg <- svyglm(smoking ~ snp + sex + colgrad +age, family=binomial(link="logit"),design = svydesign(~ 1, weights = ~ w,data=x))
summary(mreg)
?mediate
output.mediation <- mediate(mreg,yreg,treat="snp",mediator="smoking",control.value = 0,treat.value = 1)
summary(output.mediation)
