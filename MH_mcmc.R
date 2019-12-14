library(MASS)
library(rgl) # interactive 3-D plots
library(copula)
library(VineCopula)
library("scatterplot3d")
library(compositions)
library(RCurl)
library (readr)
library(repmis)
library(Hmisc)
library(Copula.Markov.survival)
library(MHadaptive)
set.seed(500)

##copula function
copulaBB8<-function(theta,delta)
{
  my_dist <- mvdc(BB8Copula(par = c(theta, delta)), 
                  margins = c("weibull","weibull"), 
                  paramMargins = list(list(shape = 2.254, scale = 6.922), 
                                      list(shape = 2.428, scale = 7.79)))
  v <- rMvdc(8570, my_dist)
  pdf_mvd <- dMvdc(v, my_dist)
  return(pdf_mvd)
}



## Define a Bayesian linear regression model
BB8_posterior<-function(pars,data)
{
  th<-pars[1] #theta
  dl<-pars[2] #delta
  
  BB8params <- prior_thndl(pars)
  likelihood_2 <- copulaBB8(BB8params[1],BB8params[2])
  posterior <- likelihood_2*prior_BB8(pars)
  return(sample(posterior, 1, replace=TRUE))
}


## Define the Prior distributions: P(theta, delta| Data)
prior_BB8<-function(pars)
{
  th<-pars[1] #theta
  dl<-pars[2] #delta

  prior_th<-runif(1,min=1,max=10)
  prior_dl<-runif(1,min=0.00001,max=0.9999)
  #prior_dl<-rgamma(1, shape=10, scale = 1)
  priors<-cbind(prior_th,prior_dl)
  
  likelihood_1 <- copulaBB8(2.27,0.9)
  prior <- likelihood_1*prior_th*prior_dl
  return(sample(prior, 1, replace=TRUE))
}


## Define the Prior distributions (theta and delta to update likelihood in posterior)
prior_thndl<-function(pars)
{
  th<-pars[1] #theta
  dl<-pars[2] #delta
  
  prior_th<-runif(1,min=1,max=10)
  prior_dl<-runif(1,min=0.00001,max=0.9999)
  #prior_dl<-rgamma(1, shape=10, scale = 1)
  priors<-cbind(prior_th,prior_dl)
  
  return(priors)
}


# simulate data
x<-rweibull(30,2.254,6.922)
y<-rweibull(30,2.428, 7.79)
u <- pobs(as.matrix(cbind(x,y)))[,1]
v <- pobs(as.matrix(cbind(x,y)))[,2]
d<-cbind(u,v)

mcmc_r<-Metro_Hastings(li_func=BB8_posterior,pars=c(1,1),
                       par_names=c('th','dl'),data=d)

## For best results, run again with the previously
## adapted variance-covariance matrix.
mcmc_r<-Metro_Hastings(li_func=li_reg,pars=c(1,1),
                       prop_sigma=mcmc_r$prop_sigma,par_names=c('th','dl'),data=d)
mcmc_r<-mcmc_thin(mcmc_r)

plotMH(mcmc_r)
