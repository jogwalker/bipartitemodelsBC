# update model

arguments <- commandArgs(T)
outdir <- arguments[1]
iter <- as.numeric(arguments[2])
rep <- as.numeric(arguments[3])
model2update <- as.character(arguments[4])


##
library(R2jags)
library(rjags)
rjags::load.module('dic')

##
load(paste(outdir,"/",model2update,sep=""))
recompile(output)
#new.output <- update(model2update,n.iter=iter)
auto <- autojags(output,n.iter=iter,n.update=rep)

save(auto,file=paste(outdir,"/",model2update,"updated.RData",sep=""))
