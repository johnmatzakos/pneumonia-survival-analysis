#Pneumonia Survival Analysis

#Author: Ioannis Matzakos Chorianopoulos

install.packages("KMsurv")
require(survival)
require(KMsurv)
data(pneumon)
#attach(pneumon)

#Generate survival estimates and plots of the Kaplan-Meier productlimit estimator for different numbers of siblings.
#creating the status variable and including it in the dataset
pneumon$status = (pneumon$chldage<12) * 1 + 0
pneumon$status

#recoding nsibs variable so that any number greater than 3 will get the value of 3
numsibs = pneumon$nsibs
numsibs
pneumon$nsibs[pneumon$nsibs>=3] = 3
pneumon$nsibs

#recoding wmonth, ages 4-6 months into value 4 and ages 7+ months into value 5
wage = pneumon$wmonth
wage
pneumon$wmonth[pneumon$wmonth>=4 & pneumon$wmonth<=6] = 4
pneumon$wmonth[wage>=7] = 5
pneumon$wmonth

#the original values of pneumon$nsibs and pneumon$wmonth are assigned to new attributes outside the dataset,
#numsibs and wage respectively, just in case they are needed later for any reason.

#fitting the model
model = survfit(Surv(chldage, status) ~ nsibs, data=pneumon, type = "kaplan")
summary(model)

#ploting the fitted model
par(mfrow=c(1,1))
plot(model, main="Kaplan-Meier Plot", xlab="Time", ylab="Probability of Survival", col=c(1,2,3,4), lwd=1)
legend('bottomright', c("nsibs=0","nsibs=1","nsibs=2","nsibs=3"), lty=1, col=c(1,2,3,4))

#log-rank test to determine whether there is a statistically significant 
#difference in the hazard rates of pneumonia given different numbers of siblings.
#H0: no difference in hazard rates between groups.
#H1: there is a difference.

#log-rank test using survdiff
survdiff(Surv(chldage, status) ~ nsibs, rho = 0, data=pneumon)

#log-rank test using the summary of a cox regression model
summary(coxph(Surv(chldage, status) ~ nsibs, data=pneumon))

#multiple comparisons using bonferroni correction
alpha = 0.05/2/6
alpha

summary(coxph(Surv(chldage, status) ~ as.factor(nsibs),data=pneumon[pneumon$nsibs<2,]))
survdiff(Surv(chldage, status) ~ as.factor(nsibs), rho = 0, data=pneumon[pneumon$nsibs<2,])
p01=5.55e-05
p01

survdiff(Surv(chldage, status) ~ as.factor(nsibs), rho = 0, data=pneumon[(pneumon$nsibs==1|pneumon$nsibs==2),])
p12=0.0826
p12

survdiff(Surv(chldage, status) ~ as.factor(nsibs), rho = 0, data=pneumon[(pneumon$nsibs==1|pneumon$nsibs==3),])
p13=0.376
p13

survdiff(Surv(chldage, status) ~ as.factor(nsibs), rho = 0, data=pneumon[(pneumon$nsibs==2|pneumon$nsibs==3),])
p23=0.842
p23

survdiff(Surv(chldage, status) ~ as.factor(nsibs), rho = 0, data=pneumon[(pneumon$nsibs==0|pneumon$nsibs==2),])
p02=2.83e-06
p02

survdiff(Surv(chldage, status) ~ as.factor(nsibs), rho = 0, data=pneumon[(pneumon$nsibs==0|pneumon$nsibs==3),])
p03=p= 0.00662
p03

sig = (p01>alpha) # FALSE
sig
sig = (p12>alpha) # TRUE
sig
sig = (p13>alpha) # TRUE
sig
sig = (p23>alpha) # TRUE
sig
sig = (p03>alpha) # FALSE
sig
sig = (p02>alpha) # FALSE
sig

#Survival analysis using Cox regression
#Backward Elimination using 0.05 significance level
#Exploratory analysis of the 2-way interactions as possible predictors for pneumonia: 
#the number of siblings, weaning age, maternal age, race, poverty, birthweight and maternal smoking.

model1 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(race)+factor(poverty)+factor(bweight)+factor(smoke))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight)+factor(smoke))+
                 mthage*(factor(race)+factor(poverty)+factor(bweight)+factor(smoke))+
                 factor(race)*(factor(poverty)+factor(bweight)+factor(smoke))+
                 factor(poverty)*(factor(bweight)+factor(smoke))+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model1)

#remove mthage:factor(poverty) 0.977833
model2 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(race)+factor(poverty)+factor(bweight)+factor(smoke))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight)+factor(smoke))+
                 mthage*(factor(race)+factor(bweight)+factor(smoke))+
                 factor(race)*(factor(poverty)+factor(bweight)+factor(smoke))+
                 factor(poverty)*(factor(bweight)+factor(smoke))+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model2)

#remove wmonth:factor(smoke) 0.907920 
model3 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(race)+factor(poverty)+factor(bweight)+factor(smoke))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                 mthage*(factor(race)+factor(bweight)+factor(smoke))+
                 factor(race)*(factor(poverty)+factor(bweight)+factor(smoke))+
                 factor(poverty)*(factor(bweight)+factor(smoke))+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model3)

#remove nsibs:factor(smoke) 0.834793
model4 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(race)+factor(poverty)+factor(bweight))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                 mthage*(factor(race)+factor(bweight)+factor(smoke))+
                 factor(race)*(factor(poverty)+factor(bweight)+factor(smoke))+
                 factor(poverty)*(factor(bweight)+factor(smoke))+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model4)

#remove nsibs:factor(race) 0.685935 
model5 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty)+factor(bweight))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                 mthage*(factor(race)+factor(bweight)+factor(smoke))+
                 factor(race)*(factor(poverty)+factor(bweight)+factor(smoke))+
                 factor(poverty)*(factor(bweight)+factor(smoke))+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model5)

#remove mthage:factor(bweight) 0.682012
model6 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty)+factor(bweight))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                 mthage*(factor(race)+factor(smoke))+
                 factor(race)*(factor(poverty)+factor(bweight)+factor(smoke))+
                 factor(poverty)*(factor(bweight)+factor(smoke))+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model6)

#remove factor(poverty):factor(race) 0.680383 
model7 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty)+factor(bweight))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                 mthage*(factor(race)+factor(smoke))+
                 factor(race)*(factor(bweight)+factor(smoke))+
                 factor(poverty)*(factor(bweight)+factor(smoke))+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model7)

#remove factor(poverty):factor(bweight) 0.576018 
model8 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty)+factor(bweight))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                 mthage*(factor(race)+factor(smoke))+
                 factor(race)*(factor(bweight)+factor(smoke))+
                 factor(poverty)*factor(smoke)+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model8)

#remove mthage:factor(race) 0.524433 
model9 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty)+factor(bweight))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                 mthage*factor(smoke)+
                 factor(race)*(factor(bweight)+factor(smoke))+
                 factor(poverty)*factor(smoke)+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model9)

#remove factor(poverty):factor(smoke) 0.474097
model10 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty)+factor(bweight))+
                 wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                 mthage*factor(smoke)+
                 factor(race)*(factor(bweight)+factor(smoke))+
                 factor(poverty)+
                 factor(bweight)*factor(smoke), data=pneumon)
anova(model10)

#remove mthage:factor(smoke) 0.464248
model11 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty)+factor(bweight))+
                  wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                  mthage+
                  factor(race)*(factor(bweight)+factor(smoke))+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model11)

#remove nsibs:factor(bweight) 0.447708
model12 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty))+
                  wmonth*(mthage+factor(race)+factor(poverty)+factor(bweight))+
                  mthage+
                  factor(race)*(factor(bweight)+factor(smoke))+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model12)

#remove wmonth:factor(race) 0.340012
model13 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty))+
                  wmonth*(mthage+factor(poverty)+factor(bweight))+
                  mthage+
                  factor(race)*(factor(bweight)+factor(smoke))+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model13)

#remove factor(race):factor(smoke) 0.241391
model14 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage+factor(poverty))+
                  wmonth*(mthage+factor(poverty)+factor(bweight))+
                  mthage+
                  factor(race)*factor(bweight)+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model14)

#remove nsibs:factor(poverty) 0.211686
model15 = coxph(Surv(chldage,status) ~ nsibs*(wmonth+mthage)+
                  wmonth*(mthage+factor(poverty)+factor(bweight))+
                  mthage+
                  factor(race)*factor(bweight)+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model15)

#remove nsibs:wmonth 0.152700 
model16 = coxph(Surv(chldage,status) ~ nsibs*mthage+
                  wmonth*(mthage+factor(poverty)+factor(bweight))+
                  mthage+
                  factor(race)*factor(bweight)+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model16)

#remove wmonth:factor(poverty) 0.099467  
model17 = coxph(Surv(chldage,status) ~ nsibs*mthage+
                  wmonth*(mthage+factor(bweight))+
                  mthage+
                  factor(race)*factor(bweight)+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model17)

#remove nsibs:mthage 0.061513 
model18 = coxph(Surv(chldage,status) ~ nsibs+
                  wmonth*(mthage+factor(bweight))+
                  mthage+
                  factor(race)*factor(bweight)+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model18)

#remove wmonth:mthage 0.053040 
model19 = coxph(Surv(chldage,status) ~ nsibs+
                  wmonth*factor(bweight)+
                  mthage+
                  factor(race)*factor(bweight)+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model19)

#remove wmonth:factor(bweight) 0.124220  
model20 = coxph(Surv(chldage,status) ~ nsibs+
                  wmonth+
                  mthage+
                  factor(race)*factor(bweight)+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model20)

#remove factor(race):factor(bweight) 0.071561  
model21 = coxph(Surv(chldage,status) ~ nsibs+
                  wmonth+
                  mthage+
                  factor(race)+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model21)

#remove factor(race) 0.3613858
model22 = coxph(Surv(chldage,status) ~ nsibs+
                  wmonth+
                  mthage+
                  factor(poverty)+
                  factor(bweight)*factor(smoke), data=pneumon)
anova(model22)

#remove factor(poverty) 0.1362200
model23 = coxph(Surv(chldage,status) ~ nsibs+wmonth+mthage+factor(bweight)*factor(smoke), data=pneumon)
anova(model23)
summary(model23)

#All the variables in the model are significant so the backward elimination process stops. 
#Additionally, the variables involved in the factor(bweight):factor(smoke) interaction can not be removed,
#because their interaction is significant. So the final model is model23

#checking the proportionality assumption
cox.zph(model23)

par(mfrow=c(2,4))
plot(cox.zph(model23), col=c(2,4))
contrasts(factor(pneumon$nsibs))
