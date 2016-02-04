
#### RUN ZIP with host species level covariate: ruminant? (for BC): Model 7

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

covars <- read.csv('~/bipartitemodelsBC/data/hostlevel.csv', head=T)
covars <- covars[with(covars,order(as.character(Host.Species))),]

str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$hostspchar <- as.numeric(covars$Ruminant)

save(long,file=paste(outdir,"/M1ilong.RData",sep=""))

## Define model

cat('      model
{    
  for (i in 1:Nobs)  #for all observations run loop around 
  { 
    count[i] ~ dpois(mu[i])  #the count data with a Poission distribution (0, 1)
    mu[i] <- lambda[i] * use[host.sp[i], par.sp[i]]  #mean abundance = characteristic abundance of each insect sp. * use of the pl.sp. k that indiv. i belongs to 
    log(lambda[i]) <- beta[par.sp[i]]      # characteristic abundance of inv.sp. j on pl.ind. i (NOTE: insect has the same abundance on all plant species)
  }   
  
  for (n in 1: Nhost.sp) #loop around all plant species, calculate use, alpha and insect diversity
  {
    for (j in 1:Npar) #loop around all insect species
    {  
      use[n, j] ~ dbern(p[n, j])  # probability of use drawn from a Bernoullli distribution 
      logit(p[n, j]) <- alpha[n, j] + alpha_cov*hostspchar[n] # p is pi, logit is link function of Bernoulli distribution
      alpha[n, j] ~ dnorm(mn[1], prec[1]) # coefficient for use is drawn from normal distribution with mean and precision)
    }
    PD_host[n]<- sum (use[n, ]) #Insect diversity = sum of the uses for each plant species
  } 
  
  # priors?   
  for (j in 1: Npar)  #loop around insect species, calculate charact abundance and host breadth
  {        
    beta[j] ~ dnorm(mn[2], prec[2]) # average abundance across all invertebrate sp. beta0 and a random effect coeff. beta1
    HB_invert[j] <- sum(use[, j]) #
  } 
  
  for (k in 1:2) # 
  { 
    mn[k] ~ dnorm(0, .0001) 
    prec[k] <- pow(sd[k], -2)
    sd[k] ~ dt(0, .016, 1)T(0, ) # a half-cauchy is apparently the t-distribution with 1 degree of freedom. 
    #sd[k] ~ dunif(0, 100)
  }
alpha_cov ~ dnorm(0, .0001)
} '    , file=(modelfile <- "~/bipartitemodelsBC/finalmodels/jags/m1i.txt"))     


## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar*long$Nhost.sp,mean=-5),ncol=long$Npar,byrow=T),
    alpha_cov=rnorm(1),
    beta=rnorm(long$Npar),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}


## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta', 'alpha','alpha_cov'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/M1ioutput.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/M1iprintoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()
