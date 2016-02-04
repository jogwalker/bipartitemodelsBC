
#### RUN ZIP with host species level covariate: feeder type (for BC): Model 8

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
#long$grazer <- as.numeric(covars$Feeder=="Grazer")
long$browser <- as.numeric(covars$Feeder=="Browser")
long$mixed <- as.numeric(covars$Feeder=="Mixed feeder")

save(long,file=paste(outdir,"/M1klong.RData",sep=""))

## Define model


modelfile <- "~/bipartitemodelsBC/finalmodels/jagsNB/m1k.txt"   


## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar*long$Nhost.sp,mean=-5),ncol=long$Npar,byrow=T),
    #alpha_g=rnorm(1),
    alpha_b=rnorm(1),
    alpha_m=rnorm(1),
    beta=rnorm(long$Npar),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}


## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta', 'alpha','alpha_m','alpha_b'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/M1koutput.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/M1kprintoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()





