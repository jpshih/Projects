#Predict number of apps using the College dataset in the ISLR package
# 1. Initial data exploration and assessment
# 2. Examine extreme values using Cook's distance
# 3. Heteroskedasticity, multicollinearity, and interaction terms
# 4. Variable selection
# 5. Dimension reduction 
# 6. Variable selection continued
# 7. Final model


#install as needed and load required packages 
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("ISLR","car","MASS","leaps","pls")
ipak(packages)

View(College)
attach(College)

# ------- 1. Initial data exploration and assessment ---------------

#initial data exploration
summary(College)
#note max. min. values that are extreme dif. from the mean (e.g., Enroll and Expend)

#build a simple linear regression and assess
reg1 = lm(Apps~.,data=College)
summary(reg1)

opar = par()
par(mfrow=c(2,2))
plot(reg1)
par(opar)

#The Residuals vs. Fitted plot shows heteroskedasticity, variance uneven
#There may be potential extreme values (e.g., Rutgers at New Brunswick)


#-------- 2. Examine extreme values ---------------

#Cook's Distance formula
cooks = function(reg,Data)
{
  n = nrow(Data)
  cutoff = 4/n #set cook's dist. cutoff
  z = round(cooks.distance(reg),4)
  View(z[z>cutoff])
  length(z[z>cutoff])
  plot(reg1, which=4,cook.levels=cutoff) #return plot
  abline(h=cutoff,col='red') 
  
}

cooks(reg1,College) 
#returns a list of colleges with extreme values
#colleges above cutoff including Rutgers NB, Harvard, and Penn State Main

#remove the rows above cutoff
College2 = College
College2$D = cooks.distance(reg1)
College2 = College2[College2$D<4/nrow(College2),]
reg2 = lm(Apps~.-D,data=College2)
z1 = round(coefficients(reg1),4)
z2 = round(coefficients(reg2),4)
print(cbind(z1,z2)) #we observe some changes in the coefficients

#An examination of the new model shows improvement
opar = par()
par(mfrow=c(2,2))
plot(reg2)
par(opar)
summary(reg2)


#Because there are extreme values, I use a robust regression (in MASS package) 
#rlm will give observations with greater Cook's distance less weight than close observations 
reg3 = rlm(Apps~.-D,data=College2)
opar = par()
par(mfrow=c(2,2))
plot(reg3)
par(opar)
summary(reg3)

#I will also try weighted least squares, which give less weight to extreme observations
z=cooks.distance(reg2)
n = nrow(College2)
cutoff = 4/length(n)
w = ifelse(z<=cutoff,1,cutoff/z)
reg4 = lm(Apps~.,weights= w,data=College2)
summary(reg4)

opar = par()
par(mfrow=c(2,2))
plot(reg4)
par(opar)
summary(reg4)


#compare the 4 models so far, reg4 with weighted least squares is the best
AIC(reg2,reg3,reg4)


# --------- 3. heteroskedasticity, multicollinearity, interaction terms ---------


#Examine heteroskedasticity using ncvTest
ncvTest(reg4)
#null hypothesis: variance is constant
#p value < 0.5, reject null hypothesis - heteroskedasticity exists
#heteroskedasticity is also observed in the "fan" shape of the Residuals vs. Fitted plot


#check if log transform is needed
reg_cox = lm(Apps~.-D,data=College2)
boxCox(reg_cox,family='yjPower',plotit=TRUE)
#lambda value close to 0, transform y

reg5 = lm(log(Apps)~.,weights= w,data=College2)
AIC(reg4,reg5) #AIC improves greatly


#should we add interaction terms?
reg_step = lm(log(Apps)~.,data=College)
reg_step = step(reg_step,~.^2)
reg_step$anova 
#Accept:Enroll is an interaction term that would greatly reduce AIC
#Accept:PhD also an interaction term, but less influential

#add these terms and examine AIC
reg6 = lm(log(Apps)~.+Accept:Enroll-D,weights= w,data=College2)
AIC(reg6,reg5) #AIC improves greatly 

reg7 = lm(log(Apps)~.+Accept:Enroll+Accept:PhD-D,weights= w,data=College2)
AIC(reg7,reg6) #AIC improves greatly after adding both terms


opar = par()
par(mfrow=c(2,2))
plot(reg7)
par(opar)
summary(reg7)


summary(reg7)
#Examine multi-collinearity
vif(reg7)
#Accept, Enroll, F.Undergrad are among some variables that exhibit VIF > 4, suggesting multicollinearity
#I will try to remove multicollinearity during the variable selection process



# ------ 4. Variable selection with leaps package ----------------------

#dataset has too many variables, we must determine which ones to focus on for the model
#explore which variables have more influence on regression accuracy

stepAIC(reg7,direction="backward") #8 variables retained
stepAIC(reg7,direction="forward") #8 variables retained

#use regsubsets to analyze which subsets have highest predictive power
regfit = regsubsets(log(Apps)~.+Accept:Enroll+Accept:PhD-D,data=College2,nvmax=15) 
summary(regfit)
#results show that variables like Accept, Accept:Enroll, Top10perc have higher predictive power

#model with 11 variables has best BIC
which.min(summary(regfit)$bic)


# ---- 5. Dimension reduction with pls package -------------


reg_pls = plsr(Apps~.+Accept:Enroll+Accept:PhD,weights=w,scale=TRUE,data=College)
summary(reg_pls)
#returns too many components

reg_pls = plsr(Apps~Accept+Top10perc-D,weights= w,data=College2)
summary(reg_pls)
#Accept and Top10Perc explain 96% of the variance
#adding more only increments the variance explanation by a little bit


# -------- 6. Variable selection cont. -----------------------------
 

#Knowing all this, let's build a model with select variables (from stepwise selection and regsubsets analysis) 
reg_8var = lm(log(Apps)~Accept+Enroll+Top10perc+Top25perc+F.Undergrad+P.Undergrad+Outstate+Accept:Enroll,weights= w,data=College2)
summary(reg_8var)
vif(reg_8var) #there is multicollinearity if using all 8 var.

#Let's try additional models, using vif to observe where multicollinearity occurs and remove it
reg_2var = lm(log(Apps)~Accept+Top10perc,weights=w,data=College2)
summary(reg_2var)
vif(reg_2var) #no multicollinearity

reg_4var = lm(log(Apps)~Accept+Top10perc+Outstate+P.Undergrad,weights=w,data=College2)
summary(reg_4var)
vif(reg_4var) #no multicollinearity

#vif analysis tells me:
#Top10Perc correlated with Top25Perc, don't need both in the model
#Accept correlated with Enroll
#F.Undergrad correlated with Accept

#compare AIC of models, reg 7 still superior
AIC(reg_8var,reg_4var,reg7)

#I will build on these models in subsequent analysis


# -------- 7. Final model: remove heteroskedascity by transforming dependent variables --------


#run boxTidwell to determine transformation values
boxTidwell(log(Apps)~Accept+Top10perc+P.Undergrad+Accept:Enroll,data=College2)

#build regression with transformed variables
reg8 = update(reg_4var,~.+I(Top10perc^3)+I(Accept^0.06)+I(P.Undergrad^0.06),data=College2)
reg9 = update(reg_2var,~.+I(Top10perc^3)+I(Accept^0.06),data=College2)
summary(reg8)
summary(reg9)

AIC(reg9,reg8,reg7) #reg8 is slightly more accurate, but it includes non-significant terms

#therefore, I select reg9 as my model


ncvTest(reg9) #p value >0.5, no heteroskedasticity
vif(reg9) #VIF close to or < 4, little to no multicollinearity

#final model:
#lm(formula = log(Apps) ~ Accept + Top10perc + I(Top10perc^3) + I(Accept^0.06), data = College2, weights = w)
#All variables significant
#with an Adjusted R-Square value of 0.97

opar = par()
par(mfrow=c(2,2))
plot(reg9)
par(opar)

summary(reg9)