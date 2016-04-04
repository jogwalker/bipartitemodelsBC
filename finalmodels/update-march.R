# 1 April 2016
# This script updates the final trunc 300k output with the aim of getting a well distributed sample of the chain
# Original run after burn in has iter 150k thinned by 150, want to discard first third or half. 
# Run update for 75k & 150k with 150 burn in - will get 500+ new samples from the chain and can discard up to 500 from old
# equivalent to running model for 375k or 450k steps with 150k+75k burn in.

## Set up output files
arguments <- commandArgs(T)
outdir <- arguments[1]
iter <- as.numeric(arguments[2])
thin <- as.numeric(arguments[3])

## load packages
library(R2jags)
library(rjags)

## load model file
load("~/bipartitemodelsBC/results/finaltrunc/cnj_output-trunc3e+05.RData")
recompile(output)
# add dic module because failed at end last time (4 April)
load.module("dic")
## update model
new.output <- update(output,n.iter=iter,n.thin=thin,progress.bar="none")

save(new.output,file=paste(outdir,"/cnj_output-trunc3e+05-update.RData")

# calculate convergence
library(jagstools)
library(dplyr)

notconv <- rhats(new.output) %>% subset(. >= 1.1) %>% length()
params <- length(rhats(new.output))


options(max.print=500000)
sink(file=paste(outdir,"/cnj_printoutput-trunc-update",iter,".txt",sep=""))
paste("not converged =", notconv, sep=" ")
paste("total params =", params, sep=" ")
print("which not converged: ")
rhats(new.output) %>% subset(. >= 1.1)
print(new.output)
sink()

