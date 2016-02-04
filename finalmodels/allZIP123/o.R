#### RUN ZIP with individual level covars (juvenile) (for BC): Model 13

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

juvenile <- as.numeric(covars$Age.class=="Juvenile")
missing <- which(is.na(juvenile)) # 35% missing

str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$juvenile <- juvenile
long$missing.ind <- missing
long$ind <- rep(1:(length(long$host.sp)/(long$Npar)),long$Npar)
long$prop <- sum(long$juvenile,na.rm=T)/length(long$juvenile)

save(long,file=paste(outdir,"/M1olong.RData",sep=""))

## Define model

cat('      model 
{    
  for (i in 1:Nobs)  #for all observations run loop around 
  { 
    count[i] ~ dpois(mu[i])  #the count data with a Poission distribution (0, 1)
    mu[i] <- lambda[i] * use[host.sp[i], par.sp[i]]  #mean abundance = characteristic abundance of each insect sp. * use of the pl.sp. k that indiv. i belongs to 
    log(lambda[i]) <- beta[par.sp[i]]  + beta_j*juvenile[ind[i]]    # characteristic abundance of inv.sp. j on pl.ind. i (NOTE: insect has the same abundance on all plant species)
  }   
  
  #Specify priors for missing data. chance of being juvenile matches prop of juvenile in rest of data
  for(i in 1:length(missing.ind)) {
    juvenile[missing.ind[i]] ~ dbern(prop)
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
  
  # priors?   
  for (j in 1: Npar)  #loop around insect species, calculate charact abundance and host breadth
  {        
    beta[j] ~ dnorm(mn[2], prec[2]) # average abundance across all invertebrate sp. beta0 and a random effect coeff. beta1
    HB_invert[j] <- sum(use[, j]) #
  } 
  
  beta_j ~ dnorm(0, 0.0001)
  
  for (k in 1:2) # 
  { 
    mn[k] ~ dnorm(0, .0001) 
    prec[k] <- pow(sd[k], -2)
    sd[k] ~ dt(0, .016, 1)T(0, ) # a half-cauchy is apparently the t-distribution with 1 degree of freedom. 
    #sd[k] ~ dunif(0, 100)
  }
} '    , file=(modelfile <- "~/bipartitemodelsBC/finalmodels/jags/m1o.txt")) 


## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar*long$Nhost.sp,mean=-5),ncol=long$Npar,byrow=T),
    beta=rnorm(long$Npar),
    beta_j=rnorm(1),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)#,
    #missing=rbinom(length(long$missing.ind),1,0.5)
  )
}


## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta', 'alpha','beta_j'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/M1ooutput.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/M1oprintoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
paste("missing data:",length(long$missing.ind)/length(long$juvenile),sep=" ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()