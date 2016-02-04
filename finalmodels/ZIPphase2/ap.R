#### RUN ZIP model with parasite genus info (Model 3)

## Set up output files
arguments <- commandArgs(T)
outdir <- arguments[1]
iter <- as.numeric(arguments[2])

## load packages
library(reshape2) #attach the reshape package for function "melt()"
library(R2jags)
library(rjags)

## read and format data
load("~/bipartitemodelsBC/data/finaldata.RData")
collec.lng$ID <- NULL

# covariates
load("~/bipartitemodelsBC/data/indlevel.RData")
female <- as.numeric(covars$Sex=="F")
missing <- which(is.na(covars$Sex)) # 68% missing... can I 

str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$par.gen <- as.factor(gsub("[_].*$","",as.character(levels(long$par.sp))))
long$Npar.gen <- length(levels(long$par.gen))
long$female <- female
long$missing.ind <- missing
long$ind <- rep(1:(length(long$host.sp)/(long$Npar)),long$Npar)

save(long,file=paste(outdir,"/aplong.RData",sep=""))

## Define model

cat('      model 
  {    
  for (i in 1:Nobs) { #for all observations run loop around 
    count[i] ~ dpois(mu[i])  #Poisson count data
    mu[i] <- lambda[i] * use[host.sp[i], par.sp[i]]  #mean abundance = abundance of each plant sp. * use of the pl.sp. k that indiv. i belongs to 
    log(lambda[i]) <- beta[par.sp[i]] + beta_s*female[ind[i]]# mean abundance of inv.sp. j on pl.ind. i
  }  
  
  for (n in 1: Nhost.sp) { #loop around all host species, calculate use, alpha and insect diversity
    
    for (j in 1:Npar) {  
      use[n, j] ~ dbern(p[n, j])  # probability of use drawn from a Bernoullli distribution 
      logit(p[n, j]) <- alpha[n, par.gen[j]] # refers to the identity list which names the category for each 
    }
    
    for (k in 1:Npar.gen) { # loop around parasite Genera
      alpha[n, k] ~ dnorm(mn[1], prec[1]) 
    }
    
    PD_host[n]<- sum (use[n, ]) #sums up the use values for each host species
  } 

    #Specify priors for missing data.
  for(i in 1:length(missing.ind)) {
    female[missing.ind[i]] ~ dbern(0.5)
  } 
  
  for (j in 1: Npar) {  #loop around parasite species, calculate charact abundance and host breadth        
    beta[j] ~ dnorm(mn[2], prec[2]) # average abundance across all invertebrate sp. beta0 and a random effect coeff. beta1
    HB_invert[j] <- sum(use[, j]) #how many hosts the parasite uses
  } 

  beta_s ~ dnorm(0, 0.0001)
  
  for (k in 1:2) # 
  { 
    mn[k] ~ dnorm(0, .0001) # mean, precision
    prec[k] <- pow(sd[k], -2)
    sd[k] ~ dt(0, .016, 1)T(0, ) # a half-cauchy is apparently the t-distribution with 1 degree of freedom. 
    #sd[k] ~ dunif(0, 100)
  }
} '    , file=(modelfile <- "~/bipartitemodelsBC/finalmodels/jags/ap.txt"))    




## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar.gen*long$Nhost.sp,mean=-5),ncol=long$Npar.gen,byrow=T),
    beta=rnorm(long$Npar),
    beta_s=rnorm(1),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}

## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta', 'alpha','beta_s'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/ap_output.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/ap_printoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()



