#### RUN MODEL with superfamily as predictor for alpha (for BC): Model 4

## Set up output files
arguments <- commandArgs(T)
outdir <- arguments[1]
iter <- as.numeric(arguments[2])


## load packages
library(R2jags)
library(rjags)

## read and format data
load("~/bipartitemodelsBC/data/finaldata.RData")
collec.lng$ID <- NULL

covars <- read.csv('~/bipartitemodelsBC/data/parasitelevel.csv', head=T)
covars <- covars[with(covars,order(Binomial)),]

str(collec.lng)
long <- as.list(collec.lng)  
long$count <- as.integer(long$count)
long$Nobs <- length(long$count)
long$Nhost.sp <- length(unique(long$host.sp))
long$Npar <- length(unique(long$par.sp))
long$par.sp <- as.factor(as.character(long$par.sp))
long$par.fam <- covars$Superfamily
long$Npar.fam <- length(unique(long$par.fam))

save(long,file=paste(outdir,"/M1blong.RData",sep=""))

## Define model


modelfile <- "~/bipartitemodelsBC/finalmodels/jagsNB/M1b.txt"        


## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar.fam*long$Nhost.sp,mean=-5),ncol=long$Npar.fam,byrow=T),
    beta=rnorm(long$Npar),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}


## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta', 'alpha'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/M1boutput.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/M1bprintoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()




