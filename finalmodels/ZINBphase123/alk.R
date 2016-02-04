#### RUN ZIP model with parasite genus info (Model 3)

## Set up output files
arguments <- commandArgs(T)
outdir <- arguments[1]
iter <- as.numeric(arguments[2])

## load packages
#library(reshape2) #attach the reshape package for function "melt()"
library(R2jags)
library(rjags)

## read and format data
load("~/bipartitemodelsBC/data/finaldata.RData")
collec.lng$ID <- NULL

covars.host <- read.csv('~/bipartitemodelsBC/data/hostlevel.csv', head=T)
covars.host <- covars.host[with(covars.host,order(as.character(Host.Species))),]

# covariates
load("~/bipartitemodelsBC/data/indlevel.RData")
rain <- covars$Precipitation
rain.sq <- sqrt(rain)
rain.norm <- (rain.sq - mean(rain.sq,na.rm=T))/sd(rain.sq,na.rm=T)
missing.rain <- which(is.na(rain))

str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$par.gen <- as.factor(gsub("[_].*$","",as.character(levels(long$par.sp))))
long$Npar.gen <- length(levels(long$par.gen))
long$rain <-  rain.norm
long$missing.ind.rain <- missing.rain
long$ind <- rep(1:(length(long$host.sp)/(long$Npar)),long$Npar)
long$browser <- as.numeric(covars$Feeder=="Browser")
long$mixed <- as.numeric(covars$Feeder=="Mixed feeder")

save(long,file=paste(outdir,"/alklong.RData",sep=""))

## Define model


modelfile <- "~/bipartitemodelsBC/finalmodels/jagsNB/alk.txt"    




## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar.gen*long$Nhost.sp,mean=-5),ncol=long$Npar.gen,byrow=T),
    beta=rnorm(long$Npar),
    beta_pr=rnorm(1), 
    alpha_b=rnorm(1),
    alpha_m=rnorm(1),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}

## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta', 'alpha','beta_pr','alpha_m','alpha_b'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/alk_output.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/alk_printoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
#paste("missing data:",length(long$missing.ind)/length(long$rain),sep=" ")
print("which not converged: ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()



