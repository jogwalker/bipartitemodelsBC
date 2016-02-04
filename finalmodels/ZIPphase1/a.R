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

str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$par.gen <- as.factor(gsub("[_].*$","",as.character(levels(long$par.sp))))
long$Npar.gen <- length(levels(long$par.gen))

save(long,file=paste(outdir,"/M1along.RData",sep=""))

## Define model

cat('      model 
  {    
  for (i in 1:Nobs) { #for all observations run loop around 
    count[i] ~ dpois(mu[i])  #Poisson count data
    mu[i] <- lambda[i] * use[host.sp[i], par.sp[i]]  #mean abundance = abundance of each plant sp. * use of the pl.sp. k that indiv. i belongs to 
    log(lambda[i]) <- beta[par.sp[i]] # mean abundance of inv.sp. j on pl.ind. i
  }  
  
  for (n in 1: Nhost.sp) { #loop around all host species, calculate use, alpha and insect diversity
    
    for (j in 1:Npar) {  
      use[n, j] ~ dbern(p[n, j])  # probability of use drawn from a Bernoullli distribution 
      logit(p[n, j]) <- alpha[n, par.gen[j]] # refers to the identity list which names the category for each 
    }
    
    for (k in 1:Npar.gen) { # loop around parasite Genera
      alpha[n, k] ~ dnorm(mn[1], prec[1]) 
    }
    
    D_byhost[n]<- sum (use[n, ]) #sums up the use values for each host species
  } 
  
  for (j in 1: Npar) {  #loop around parasite species, calculate charact abundance and host breadth        
    beta[j] ~ dnorm(mn[2], prec[2]) # average abundance across all invertebrate sp. beta0 and a random effect coeff. beta1
    HB_invert[j] <- sum(use[, j]) #how many hosts the parasite uses
  } 
  
  for (k in 1:2) # 
  { 
    mn[k] ~ dnorm(0, .0001) # mean, precision
    prec[k] <- pow(sd[k], -2)
    sd[k] ~ dt(0, .016, 1)T(0, ) # a half-cauchy is apparently the t-distribution with 1 degree of freedom. 
    #sd[k] ~ dunif(0, 100)
  }
} '    , file=(modelfile <- "~/bipartitemodelsBC/finalmodels/jags/M1a.txt"))    




## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar.gen*long$Nhost.sp,mean=-5),ncol=long$Npar.gen,byrow=T),
    beta=rnorm(long$Npar),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}

## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'D_byhost', 'beta', 'alpha'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/M1aoutput.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/M1aprintoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()



