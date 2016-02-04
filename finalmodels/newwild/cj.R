
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

covars.host <- read.csv('~/bipartitemodelsBC/data/hostlevel-switchwild.csv', head=T)
covars.host <- covars.host[with(covars.host,order(as.character(Host.Species))),]

str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$domestic <- as.numeric(covars.host$Wild)


save(long,file=paste(outdir,"/cjlong.RData",sep=""))

## Define model


modelfile <- "~/bipartitemodelsBC/finalmodels/jagsNB/phase5/cj.txt"     


## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar*long$Nhost.sp,mean=-5),ncol=long$Npar,byrow=T),
    alpha_d=rnorm(1),
    beta=matrix(rnorm(long$Npar*long$Nhost.sp,mean=-5),ncol=long$Npar,byrow=T),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}


## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta', 'alpha','alpha_d','prec.beta','r'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/cj_output.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/cj_printoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()
