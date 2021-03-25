#===================================
# Problem Set A2 - Gregory Campbell
#===================================

# Preliminaries 

setwd("C:/Users/Gregory Campbell/Downloads/Econ 613")

rm(list = ls())

library(bayesm)
library(tidyverse)
library(mlogit)

data("margarine")

install.packages("nloptr")

library(nloptr)

#==========================================
# Exercise 1: Data Description
#==========================================

set.seed(123)

# Finding average of product characteristics

colMeans(margarine$choicePrice[3:12])

# Now, the standard deviations 

margarine$choicePrice[3:12] %>% apply(2,sd) 

# Market share 

100*table(margarine$choicePrice[2])/nrow(margarine$choicePrice)

# Market share by product characteristics 

# Finding the mean price of all products 

# Creating matrix of product prices

prices = margarine$choicePrice[3:12]

p_hat = sum(colSums(prices))/(10*nrow(prices))

# Making a single list of the products chosen and their respective prices

choice_1 = margarine$choicePrice %>% filter(choice == 1) %>% select(2,3) %>% rename(Price = PPk_Stk)
choice_2 = margarine$choicePrice %>% filter(choice == 2) %>% select(2,4) %>% rename(Price = PBB_Stk)
choice_3 = margarine$choicePrice %>% filter(choice == 3) %>% select(2,5) %>% rename(Price = PFl_Stk)
choice_4 = margarine$choicePrice %>% filter(choice == 4) %>% select(2,6) %>% rename(Price = PHse_Stk)
choice_5 = margarine$choicePrice %>% filter(choice == 5) %>% select(2,7) %>% rename(Price = PGen_Stk)
choice_6 = margarine$choicePrice %>% filter(choice == 6) %>% select(2,8) %>% rename(Price = PImp_Stk)
choice_7 = margarine$choicePrice %>% filter(choice == 7) %>% select(2,9) %>% rename(Price = PSS_Tub)
choice_8 = margarine$choicePrice %>% filter(choice == 8) %>% select(2,10) %>% rename(Price = PPk_Tub)
choice_9 = margarine$choicePrice %>% filter(choice == 9) %>% select(2,11) %>% rename(Price = PFl_Tub)
choice_10 = margarine$choicePrice %>% filter(choice == 10) %>% select(2,12) %>% rename(Price = PHse_Tub)

# The list itself

choice_list = rbind(choice_1,choice_2,choice_3,choice_4,choice_5,choice_6,choice_7,choice_8,
                   choice_9,choice_10)

# Looking at the products that were chosen whose prices are below the mean price of all products

below_avg = choice_list %>% filter(Price < p_hat)

100*table(below_avg$choice)/nrow(below_avg)

# Now, looking at products chosen whose prices are above the mean of all products

above_avg = choice_list %>% filter(Price >= p_hat)

100*table(above_avg$choice)/nrow(above_avg)


# Mapping between observed attributes and choices 

# Merging the Choice price and Demo datasets for margarine

marg_merge = left_join(margarine$choicePrice,margarine$demos)

# Looking at income bins (22.5,47.5,130)

hhincome_low = marg_merge %>% filter(Income <= 22.5)

# Low income households

100*table(hhincome_low$choice)/nrow(hhincome_low)

# So, it looks like low income households tend to choose product 1 (PPk_Stk)

# Medium income households 

hhincome_med = marg_merge %>% filter(Income > 22.5 & Income <= 47.5)

100*table(hhincome_med$choice)/nrow(hhincome_med)

# So, it looks like medium income households still prefer product 1

# High income households

hhincome_high = marg_merge %>% filter(Income > 47.5)

100*table(hhincome_high$choice)/nrow(hhincome_high)

# So, it still looks like high income households prefer product 1, but products 8 and 9 seem
# to be more popular in the high income households than the medium and low income households

# College educated

college_pref = marg_merge %>% filter(college == 1)

100*table(college_pref$choice)/nrow(college_pref)

# Looks like the college educated prefer product 1

# No college

nocollege_pref = marg_merge %>% filter(college == 0)

100*table(nocollege_pref$choice)/nrow(nocollege_pref)

# Looks incredibly similar to what college graduates buy

# Retired

retired_pref = marg_merge %>% filter(retired == 1)

100*table(retired_pref$choice)/nrow(retired_pref)

# Again, it looks like retired people prefer product 1, but also seem to like product 3 more
# than most groups

# Family size less than 3

Fsl3_pref = marg_merge %>% filter(Fam_Size < 3)

100*table(Fsl3_pref$choice)/nrow(Fsl3_pref)

# It looks like smaller families prefer product 1

# Family size of 3 or 4

Fs3_4_pref = marg_merge %>% filter(Fs3_4 == 1)

100*table(Fs3_4_pref$choice)/nrow(Fs3_4_pref)

# Again, the main choices of families of 3 to 4 is product 1, but it also like these families 
# aren't fond of product 9 compared to the smaller sized households

# Family size 5 or above

Fs5_pref = marg_merge %>% filter(Fs5. == 1)

100*table(Fs5_pref$choice)/nrow(Fs5_pref)

# Again, it looks like large families prefer product 1 and also like producy 4


#==========================================
# Exercise 2: First Model
#==========================================

# For this model, I'm going to use a conditional logit model to see how prices affect which
# of the 10 products are selected by individuals. The reason I'm using a conditional logit is
# because we're trying to measure the effect of price on demand for a product and price is 
# a product level characteristic. Therefore, a conditional model makes the most sense and something
# that I am assuming is that the changes in price of a product will have roughly the same effect
# across all products, which seems reasonable to me given that price tends to have a negative effect
# on demand in general. So, we are going to regress the product choices people make on the prices of
# these products. 


choice = margarine$choicePrice$choice

like_cond = function(param,choice,prices)
  {
    ni = nrow(prices)
    nj = ncol(prices)
    ut = mat.or.vec(ni,nj)
    
    ut[,1] = param[10]*prices[,1]

    for (j in 2:nj)
    {
      # conditional logit
      
      ut[,j] = param[j-1] + param[10]*prices[,j]
  
    }
    prob   = exp(ut)            
    prob   = sweep(prob,MARGIN=1,FUN="/",STATS=rowSums(prob)) 
    
    probc = NULL
    
    for (i in 1:ni)
    {
      probc[i] = prob[i,choice[i]]
    }
    probc[probc>0.999999] = 0.999999
    probc[probc<0.000001] = 0.000001
    like = sum(log(probc))
    return(-like)
  }

npar = 10 
lower  = rep(-10,npar)
upper  = rep(10,npar)
start  = runif(npar) 

res_logit_sub   = nloptr(start,eval_f=like_cond, lb=lower,ub=upper,
                          opts=list("algorithm"="NLOPT_LN_SBPLX","xtol_rel"=1.0e-10, print_level = 0,"maxeval"=10000),
                          choice=choice,prices=prices)

res_logit_sub$solution

res_logit_sub$solution[10]

# So, what this price variable is indicating is that an increase in the price of any of the
# 10 products decreases the likelihood of people selecting that product, as one would expect. 


#==========================================
# Exercise 3: Second Model
#==========================================

# For this problem I am going to use a multinomial logit model. The reason for this is because we want to 
# estimate the effect family income has on which products are chosen. As we saw in exercise 1, varying 
# family income has an effect on the probability of a family choosing certain products. Therefore, we
# want to track the effect income has on the probability of selecting each product because the effect
# will be different for each product. Thus, a multinomial logit is a sensible model to use since it
# is able to track these differing effects income has on product selection. Also, in this model we're
# going to control for effects family size, individuals being whitecollar, retired, or college educated
# might have on choice to get better beta estimates for income.  

Income = marg_merge$Income

Whitecollar = marg_merge$whtcollar

Retired = marg_merge$retired

Fam_s = marg_merge$Fam_Size

College = marg_merge$college

multi_fun = function(param,choice,Retired,Fam_s,Income,Whitecollar,College)
{
  ni = length(Income)
  nj = 10
  parameters = 54
  ut = mat.or.vec(ni,nj)
  
  # multinomial logit
  
  ut[,1] = 0  # Setting reference category 
  
  
  for (j in 2:nj)
  {
    # multinomial logit
    
    ut[,j] = param[j-1] + param[j+8]*Income + param[j+17]*Whitecollar + 
      param[j+26]*College + param[j+35]*Retired + param[j+44]*Fam_s
  }
  
  prob   = exp(ut)            
  prob   = sweep(prob,MARGIN=1,FUN="/",STATS=rowSums(prob)) 
  
  probc = NULL
  for (i in 1:ni)
  {
    probc[i] = prob[i,choice[i]]
  }
  probc[probc>0.999999] = 0.999999
  probc[probc<0.000001] = 0.000001
  like = sum(log(probc))
  return(-like)
}

npar= 54
lower  = rep(-10,npar)
upper  = rep(10,npar)
start  = runif(npar)

res_mult_logit   = optim(start,fn=multi_fun,method = "BFGS",control = list(trace=6,REPORT=1,maxit=10000),
                         choice = choice, Income=Income, Whitecollar=Whitecollar,College=College,
                         Fam_s=Fam_s,Retired=Retired,hessian = FALSE)

res_mult_logit$par[10:18]

# What these betas in front of family income are telling us is effect an increase in family income has on
# the likelihood of families to choose one product compared to the reference category, which is product 1
# in my situation. For example, the beta coefficient in found in front of product 2 is roughly 
# -0.00175. This means that families whose incomes increase are less likely to choose product 1 compared
# to product 2, holding all else equal. 

#==========================================
# Exercise 4: Marginal Effects
#==========================================

# Conditional logit: 

clogit_res = res_logit_sub$solution


# Computing probability matrix

ut=mat.or.vec(nrow(prices),ncol(prices))
ut[,1] = clogit_res[10]*prices[,1]
for (j in 2:ncol(prices))
{
  # conditional logit
  ut[,j] = clogit_res[j-1] + clogit_res[10]*prices[,j]
}
prob   = exp(ut)          
prob   = sweep(prob,MARGIN=1,FUN="/",STATS=rowSums(prob)) 

# Computing marginal effects

mgl_effects_1=mat.or.vec(nrow(prices),ncol(prices))

dummy_vars = c(1,rep(0,9))

for(j in 1:10){
  mgl_effects_1[,j]=prob[,j]*(dummy_vars[j]-prob[,1])*clogit_res[10]
}

avg_mgl_effects=colMeans(mgl_effects_1)

avg_mgl_effects

# So, as we can see from these marginal effects, increasing the price of product 1 decreases the likelihood
# of families choosing that product, but also increases the likelihood of families choosing other products.


# Average marginal effects for the multinomial logit. 

# Coefficients

coef_mult = res_mult_logit$par


# Making the probability matrix: 

ut=mat.or.vec(length(Income),ncol(prices))

ut[,1] = 0

for (j in 2:10){
  # multinomial logit
  
  ut[,j] = coef_mult[j-1] + coef_mult[j+8]*Income + coef_mult[j+17]*Whitecollar + 
    coef_mult[j+26]*College + coef_mult[j+35]*Retired + coef_mult[j+44]*Fam_s
}

prob   = exp(ut)

prob   = sweep(prob,MARGIN=1,FUN="/",STATS=rowSums(prob)) 

# Computing marginal effects:

# Income 

mgl_effects_Income = mat.or.vec(nrow(prices),10)

Income_coefs = c(0,coef_mult[10:18])

wgt_avg_Income = prob%*%Income_coefs

for(j in 1:10){
  mgl_effects_Income[,j]=prob[,j]*(Income_coefs[j]-wgt_avg_Income)
}

avg_mgl_effects_Income = colMeans(mgl_effects_Income)

avg_mgl_effects_Income

# So, what these marginal effects are telling us is the effect an increase in family income has on the 
# likelihood of choosing each product, holding all else constant. For example, the marginal effect of an 
# increase in family income on product 1 is roughly -0.0014. What this means is that an increase in 
# family income decreases the likelihood of households choosing product 1 by roughly 0.0014 log points. 

# Whitecollar marginal effects

mgl_effects_Whitecollar = mat.or.vec(nrow(prices),10)

Whitecollar_coefs = c(0,coef_mult[19:27])

wgt_avg_Whitecollar = prob%*%Whitecollar_coefs

for(j in 1:10){
  mgl_effects_Whitecollar[,j]=prob[,j]*(Whitecollar_coefs[j]-wgt_avg_Whitecollar)
}

avg_mgl_effects_Whitecollar = colMeans(mgl_effects_Whitecollar)

avg_mgl_effects_Whitecollar

# What these marginal effects are telling us is the effect a household having white collar workers has 
# on their likelihood of choosing each product. 

# College marginal effects

mgl_effects_College = mat.or.vec(nrow(prices),10)

College_coefs = c(0,coef_mult[28:36])

wgt_avg_College = prob%*%College_coefs

for(j in 1:10){
  mgl_effects_College[,j]=prob[,j]*(College_coefs[j]-wgt_avg_College)
}

avg_mgl_effects_College = colMeans(mgl_effects_College)

avg_mgl_effects_College

# What these marginal effects are telling us is the effect a household having college educated residents has 
# on their likelihood of choosing each product.

# Retired marginal effects

mgl_effects_Retired = mat.or.vec(nrow(prices),10)

Retired_coefs = c(0,coef_mult[37:45])

wgt_avg_Retired = prob%*%Retired_coefs

for(j in 1:10){
  mgl_effects_Retired[,j]=prob[,j]*(Retired_coefs[j]-wgt_avg_Retired)
}

avg_mgl_effects_Retired = colMeans(mgl_effects_Retired)

avg_mgl_effects_Retired

# What these marginal effects are telling us is the effect a household having retired workers has 
# on their likelihood of choosing each product.

# Family size marginal effects

mgl_effects_Fam = mat.or.vec(nrow(prices),10)

Fam_coefs = c(0,coef_mult[46:54])

wgt_avg_Fam = prob%*%Fam_coefs

for(j in 1:10){
  mgl_effects_Fam[,j]=prob[,j]*(Fam_coefs[j]-wgt_avg_Fam)
}

avg_mgl_effects_Fam = colMeans(mgl_effects_Fam)

avg_mgl_effects_Fam

# What these marginal effects are telling us is the effect an increase in family size has on their 
# likelihood of choosing each product. 

#==========================================
# Exercise 5: IIA
#==========================================

# Mixed logit:

mixed_fun = function(param,choice,prices,Retired,Fam_s,Income,Whitecollar,College)
{
  ni = length(Income)
  nj = 10
  parameters = 47
  ut = mat.or.vec(ni,nj)
  
  ut[,1] = 0  # Setting reference category 
  
  
  for (j in 2:nj)
  {
    # mixed logit
    
    ut[,j] = param[j-1] + param[10]*prices[,j] + param[11]*College + param[j+10]*Whitecollar + 
      param[j+19]*College + param[j+28]*Retired + param[j+36]*Fam_s
  }
  
  prob   = exp(ut)            
  prob   = sweep(prob,MARGIN=1,FUN="/",STATS=rowSums(prob)) 
  
  probc = NULL
  for (i in 1:ni)
  {
    probc[i] = prob[i,choice[i]]
  }
  probc[probc>0.999999] = 0.999999
  probc[probc<0.000001] = 0.000001
  like = sum(log(probc))
  return(-like)
}

npar= 47
lower  = rep(-10,npar)
upper  = rep(10,npar)
start  = runif(npar)

res_mixed_logit   = optim(start,fn=mixed_fun,method = "BFGS",control = list(trace=6,REPORT=1,maxit=10000),
                         choice = choice, Income=Income, Whitecollar=Whitecollar,College=College,
                         prices=prices,Fam_s=Fam_s,Retired=Retired,hessian = FALSE)

Beta_f = res_mixed_logit$par


# Restricted model:

# Making a new restricted dataset with household characteristics 

marg_merg2  = left_join(margarine$choicePrice,margarine$demos) %>% filter(choice < 10) %>% 
  select(-12) %>% select(-1)
  
# Making characteristic variables (you can ignore this)

choice_2 = marg_merg2$choice
prices_2 = marg_merg2[,2:10]
Income_2 = marg_merg2$Income
Fams_2 = marg_merg2$Fam_Size
College_2 = marg_merg2$college
Whitecollar_2 = marg_merg2$whtcollar
Retired_2 = marg_merg2$retired


# Estimating restricted model
  
  mixed_fun_2 = function(param,choice_2,prices_2,Income_2,Fams_2,College_2,Whitecollar_2,Retired_2)
  {
    ni = length(Income_2)
    nj = 9
    parameters = 42
    ut = mat.or.vec(ni,nj)
    
    ut[,1] = 0  # Setting reference category 
    
    
    for (j in 2:nj)
    {
      # mixed logit
      
      ut[,j] = param[j-1] + param[9]*prices_2[,j] + param[10]*College_2 + param[j+9]*Whitecollar_2 + 
        param[j+17]*Income_2 + param[j+25]*Retired_2 + param[j+33]*Fams_2
    }
    
    prob   = exp(ut)            
    prob   = sweep(prob,MARGIN=1,FUN="/",STATS=rowSums(prob)) 
    
    probc = NULL
    for (i in 1:ni)
    {
      probc[i] = prob[i,choice_2[i]]
    }
    probc[probc>0.999999] = 0.999999
    probc[probc<0.000001] = 0.000001
    like = sum(log(probc))
    return(-like)
  }
  

npar= 42
lower  = rep(-10,npar)
upper  = rep(10,npar)
start  = runif(npar)

res_mixed_logit_2   = optim(start,fn=mixed_fun_2,method = "BFGS",control = list(trace=6,REPORT=1,maxit=10000),
                          choice_2=choice_2,Income_2=Income_2, Whitecollar_2=Whitecollar_2,
                          College_2=College_2,prices_2=prices_2,Fams_2=Fams_2,
                          Retired_2=Retired_2,hessian = FALSE)


# Coefficient estimates for restricted model

Beta_r = res_mixed_logit_2$par

Beta_r 

r_likelihood = res_mixed_logit_2$value

# Compared to the previous mixed logit I ran, the beta estimates in front of my variables are clearly 
# different. Clearly, IIA is violated, which is to no ones surprise.  


# IIA test statistic:

param = Beta_f[-c(9,20,29,38,47)]

Lr_bf_like = mixed_fun_2(param,choice_2,prices_2,Income_2,Fams_2,College_2,Whitecollar_2,Retired_2)
            
mtt = -2*(Lr_bf_like - r_likelihood)

# Calculating mtt

mtt

# Based on this mtt value, it looks like we can reject the hypothesis that IIA is not violated. 

