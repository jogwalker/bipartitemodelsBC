#### RUN ZIP MODEL with host genus random effect on use (for BC): Model 18

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

set.seed(3)

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
long$hostgen <- as.factor(as.character(covars$Genus))
long$ngen <- length(levels(long$hostgen))

save(long,file=paste(outdir,"/M1glong",iter,".RData",sep=""))

## Define model


modelfile <- "~/bipartitemodelsBC/finalmodels/jagsNB/beta/m1g.txt" 


## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5,0.5), 
    sd=c(3,1.5,3),
    alpha=matrix(rnorm(long$Npar*long$Nhost.sp,mean=-5),ncol=long$Npar,byrow=T),
    alpha_gen=rnorm(long$ngen),
    beta=rnorm(long$Npar),
    use=matrix(rep(0.9,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}


## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'phi','beta','alpha','r','HB_invert','PD_host','hosts'), modelfile, n.chains=3, n.iter=iter,n.burnin=burnin,n.thin=thin) # or use defaults

save(output, file = paste(outdir,"/M1goutput",iter,".RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/M1gprintoutput",iter,".txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()
