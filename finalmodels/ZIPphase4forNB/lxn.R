#### RUN ZIP model with individual level covars (for BC): Model 10 precipitation

## Set up output files
arguments <- commandArgs(T)
outdir <- arguments[1]
iter <- as.numeric(arguments[2])

## load packages
#library(reshape2) #attach the reshape package for function "melt()"
library(R2jags)
library(rjags)


load("~/bipartitemodelsBC/data/finaldata.RData")
collec.lng$ID <- NULL

# covariates
load("~/bipartitemodelsBC/data/indlevel.RData")
rain <- covars$Precipitation
rain.sq <- sqrt(rain)
rain.norm <- (rain.sq - mean(rain.sq,na.rm=T))/sd(rain.sq,na.rm=T)
missing.rain <- which(is.na(rain))

load("~/bipartitemodelsBC/data/indlevel.RData")
treated <- as.numeric(covars$Previous.Treatment!="None" & !is.na(covars$Previous.Treatment))
missing.treat <- which(is.na(covars$Previous.Treatment))


str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$rain <-  rain.norm
long$missing.ind.rain <- missing.rain
long$treated <- treated
long$missing.ind.treat <- missing.treat
long$ind <- rep(1:(length(long$host.sp)/(long$Npar)),long$Npar)

save(long,file=paste(outdir,"/lxnlong.RData",sep=""))

## Define model

cat('      model 
{    
  for (i in 1:Nobs)  #for all observations run loop around 
  { 
    count[i] ~ dpois(mu[i])  #the count data with a Poission distribution (0, 1)
    mu[i] <- lambda[i] * use[host.sp[i], par.sp[i]]  #mean abundance = characteristic abundance of each insect sp. * use of the pl.sp. k that indiv. i belongs to 
    log(lambda[i]) <- beta[par.sp[i]]   + beta_pr[1]*rain[ind[i]] + beta_pr[2]*rain[ind[i]]*treated[ind[i]]+ beta_t*treated[ind[i]]   # characteristic abundance of inv.sp. j on pl.ind. i (NOTE: insect has the same abundance on all plant species)
  }   
  
  for (n in 1: Nhost.sp) #loop around all plant species, calculate use, alpha and insect diversity
  {
    for (j in 1:Npar) #loop around all insect species
    {  
      use[n, j] ~ dbern(p[n, j])  # probability of use drawn from a Bernoullli distribution 
      logit(p[n, j]) <- alpha[n, j] # p is pi, logit is link function of Bernoulli distribution
      alpha[n, j] ~ dnorm(mn[1], prec[1]) # coefficient for use is drawn from normal distribution with mean and precision)
    }
    
    PD_host[n]<- sum (use[n, ]) #Insect diversity = sum of the uses for each plant species
  } 
  
  #Specify priors for missing data.
  for(i in 1:length(missing.ind.rain)) {
    rain[missing.ind.rain[i]] ~ dnorm(0,1)T(-1.738577, )
  }  

  for(i in 1:length(missing.ind.treat)) {
    treated[missing.ind.treat[i]] ~ dbern(0.5)
  } 
  
  # priors?   
  for (j in 1: Npar)  #loop around insect species, calculate charact abundance and host breadth
  {        
    beta[j] ~ dnorm(mn[2], prec[2]) # average abundance across all invertebrate sp. beta0 and a random effect coeff. beta1
    HB_invert[j] <- sum(use[, j]) #
  } 
  
  beta_t ~ dnorm(0, 0.0001)
  
  for (k in 1:2) # 
  { 
    beta_pr[k] ~ dnorm(0, 0.0001)
    mn[k] ~ dnorm(0, .0001) 
    prec[k] <- pow(sd[k], -2)
    sd[k] ~ dunif(0, 100)
  }
} '    , file=(modelfile <- "~/bipartitemodelsBC/finalmodels/jags/phase4/lxn.txt"))   


## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar*long$Nhost.sp,mean=-5),ncol=long$Npar,byrow=T),
    beta=rnorm(long$Npar),
    beta_pr=rnorm(2),
    beta_t=rnorm(1),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}


## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta','beta_pr','alpha','beta_t'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/lxn_output.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/lxn_printoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
#paste("missing data:",length(long$missing.ind)/length(long$rain),sep=" ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()