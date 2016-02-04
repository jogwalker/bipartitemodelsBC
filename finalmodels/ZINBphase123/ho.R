#### RUN ZIP MODEL with host genus random effect on use (for BC): Model 19

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

covars.host <- read.csv('~/bipartitemodelsBC/data/hostlevel.csv', head=T)
covars.host <- covars.host[with(covars.host,order(as.character(Host.Species))),]

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
long$hostgen <- as.factor(as.character(covars.host$Genus))
long$ngen <- length(levels(long$hostgen))
long$juvenile <- juvenile
long$missing.ind <- missing
long$ind <- rep(1:(length(long$host.sp)/(long$Npar)),long$Npar)
long$prop <- sum(long$juvenile,na.rm=T)/length(long$juvenile)

save(long,file=paste(outdir,"/holong.RData",sep=""))

## Define model


modelfile <- "~/bipartitemodelsBC/finalmodels/jagsNB/ho.txt"    

## Define initial values

inits <- function() {
  list(
    mn=c(0.5,-5), 
    sd=c(3,1.5),
    alpha=matrix(rnorm(long$Npar*long$ngen,mean=-5),ncol=long$Npar,byrow=T),
    beta=rnorm(long$Npar),
    beta_j=rnorm(1),
    use=matrix(rep(1,long$Npar*long$Nhost.sp),ncol=long$Npar,byrow=T)
  )
}


## Run model
output  <- jags(long, inits = inits, c('mn', 'sd', 'use', 'HB_invert', 'PD_host', 'beta', 'alpha','beta_j'), modelfile, n.chains=3, n.iter=iter) # or use defaults

save(output, file = paste(outdir,"/ho_output.RData",sep=""))


# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(output))


options(max.print=100000)
sink(file=paste(outdir,"/ho_printoutput.txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
rhats(output) %>% subset(. >= 1.1)
print(output)
sink()
