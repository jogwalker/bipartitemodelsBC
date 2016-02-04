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
missing <- which(is.na(rain))

covars.host <- read.csv('~/bipartitemodelsBC/data/hostlevel.csv', head=T)
covars.host <- covars.host[with(covars.host,order(as.character(Host.Species))),]

str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$rain <-  rain.norm
long$missing.ind <- missing
long$ind <- rep(1:(length(long$host.sp)/(long$Npar)),long$Npar)
long$browser <- as.numeric(covars.host$Feeder=="Browser")
long$mixed <- as.numeric(covars.host$Feeder=="Mixed feeder")

save(long,file=paste(outdir,"/lklong.RData",sep=""))

## Define model


modelfile <- "~/bipartitemodelsBC/finalmodels/jagsNB/lk.txt"   


## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar*long$Nhost.sp,mean=-5),ncol=long$Npar,byrow=T),
    beta=rnorm(long$Npar),
    beta_pr=rnorm(1),
    alpha_b=rnorm(1),
    alpha_m=rnorm(1),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}


## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta','beta_pr','alpha','alpha_b','alpha_m'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/lk_output.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/lk_printoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
#paste("missing data:",length(long$missing.ind)/length(long$rain),sep=" ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()