#### RUN ZIP with individual level covars (for BC): Model 11

## Set up output files
arguments <- commandArgs(T)
outdir <- arguments[1]
iter <- as.numeric(arguments[2])
burnin <- as.numeric(arguments[3])
thin <- as.numeric(arguments[4])

## load packages
#library(reshape2) #attach the reshape package for function "melt()"
library(R2jags)
library(rjags)


load("~/bipartitemodelsBC/data/finaldata.RData")
collec.lng$ID <- NULL

# covariates
load("~/bipartitemodelsBC/data/indlevel.RData")
# covars <- read.csv('~/bipartitemodelsBC/data/indlevel.csv', head=T)#, stringsAsFactors=FALSE)
# # already sorted, need to repeat each within model
female <- as.numeric(covars$Sex=="F")
missing <- which(is.na(covars$Sex)) # 68% missing... can I actually use this?

str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$female <- female
long$missing.ind <- missing
long$ind <- rep(1:(length(long$host.sp)/(long$Npar)),long$Npar)

save(long,file=paste(outdir,"/M1plong",iter,".RData",sep=""))

## Define model


modelfile <- "~/bipartitemodelsBC/finalmodels/jagsNB/beta/m1p.txt"   


## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar*long$Nhost.sp,mean=-5),ncol=long$Npar,byrow=T),
    beta=rnorm(long$Npar),
    beta_s=rnorm(1),
    use=matrix(rep(0.9,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)#,
    #missing=rbinom(length(long$missing.ind),1,0.5)
  )
}


### Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'phi','mu','beta_s'), modelfile, n.chains=3, n.iter=iter,n.burnin=burnin,n.thin=thin) # or use defaults

save(output, file = paste(outdir,"/M1poutput",iter,".RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/M1pprintoutput",iter,".txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
paste("missing data:",length(long$missing.ind)/length(long$female),sep=" ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()