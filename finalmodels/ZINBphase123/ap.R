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


modelfile <- "~/bipartitemodelsBC/finalmodels/jagsNB/ap.txt"    




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



