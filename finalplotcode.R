# final plot code ch3
# update for new laptop 29 April 2016

#load data from updated run
#load("~/dat/bipartitemodelsBC/cnj_output-trunc3e+05-update.RData")
load("~/dat/bipartitemodelsBC/march-cnj_output-trunc3e+05.RData")

library(coda,quietly=TRUE,verbose=FALSE)
library(R2jags,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)

# pdf("~/dat/bipartitemodelsBC/figs/trace-update.pdf",height=10)
# traceplot(new.output,mfrow=c(4,2),varname=c('mn', 'sd','alpha_d','beta_t','r','deviance'),ask=FALSE)
# dev.off()

pdf("~/dat/bipartitemodelsBC/figs/trace-300k-long.pdf",height=10)
traceplot(output,mfrow=c(4,2),varname=c('mn', 'sd','alpha_d','beta_t','r','deviance'),ask=FALSE)
dev.off()


attach.jags(output)
pdf("~/dat/bipartitemodelsBC/figs/density.pdf",height=10)
par(mfrow = c(4,2))
plot(density(mn[,1]),main="mn[1]")
plot(density(mn[,2]),main="mn[2]")
plot(density(sd[,1]),main="sd[1]")
plot(density(sd[,2]),main="sd[2]")
plot(density(alpha_d),main="alpha_d")
plot(density(beta_t),main="beta_t")
plot(density(r),main="r")
plot(density(deviance),main="deviance")
dev.off()


load("~/git/bipartitemodelsBC/data/finaldata.RData")
load("~/dat/bipartitemodelsBC/march-cnjlong-trunc3e+05.RData")

library(dplyr,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(tidyr,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(reshape2,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)


# observed prevalence
agg <- collec.lng %>% group_by(host.sp,par.sp) %>% summarize(nhost=length(unique(ID)),totalcount=sum(count),npos=sum(as.numeric(count>0)))
agg$prev <- agg$npos/agg$nhost

# predicted probability of use
PMat <- t(apply(output$BUGSoutput$sims.list$use, c(2, 3), 
                function(x) ifelse(median(x) == 1 | median(x) == 0, abs(mean(x) - .00001), median(x)))) # summarize across the chain



dimnames(PMat) <- list(levels(long$par.sp),levels(long$host.sp))
prob.use1<- melt(PMat)
prob.use1<-prob.use1[,c(2,1,3)]
names(prob.use1)<-c( 'host.sp', 'par.sp', 'use')
prob.use<-merge(x=prob.use1, y=agg, by = c('host.sp', 'par.sp') , all=T)
prob.use[is.na(prob.use)] <- 0 

# hist(prob.use$use, main = "Estimated Probability of Use", xlab="Probability of Use")
# hist(prob.use$prev)
# library(boot)
# hist(use[,16,60])
# hist(alpha[,16,60])
# hist(inv.logit((alpha[,16,60])))
plot(prob.use$use ~ prob.use$prev)
library(ggplot2)
pdf("~/dat/bipartitemodelsBC/figs/probuse-rotate.pdf",height=4)
ggplot(data=prob.use,aes(y=prev,x=use)) + geom_point() + theme_bw() + ylab("Mean estimate for use") + xlab("Observed prevalence") 
dev.off()

print("Correlation between predicted probability of use and observed prevalence:")
cor(prob.use$use, prob.use$prev) 


# counts vs abundance
# beta is the expected abundance of each host*parasite, without treatment
CMat <- t(apply(output$BUGSoutput$sims.list$beta, c(2, 3), median))
eCMat <- exp(CMat)
dimnames(CMat) <- list(levels(long$par.sp),levels(long$host.sp))



agg$meancount <- agg$totalcount/agg$nhost
# observed counts
expabund <- melt(CMat)
expabund<-expabund[,c(2,1,3)]
names(expabund)<-c( 'host.sp', 'par.sp', 'beta')
abund<-merge(x=expabund, y=agg, by = c('host.sp', 'par.sp') , all=T)

cor(abund$beta,abund$meancount)
no0 <- abund %>% filter(meancount > 0)
cor(exp(no0$beta),no0$meancount)

plot(abund$beta ~ log(abund$meancount+0.5))

abund$logcount <- ifelse(abund$meancount==0,log(abund$meancount+1),log(abund$meancount))
abund$log <- ifelse(abund$meancount==0,TRUE,FALSE)

pdf("~/dat/bipartitemodelsBC/figs/abund-rotate.pdf",height=4)
ggplot(data=abund,aes(y=logcount,x=beta,colour=log)) + theme_bw() + ylab("Log predicted mean abundance") +xlab("Log observed mean count") + geom_point() + theme(legend.position="none")
dev.off()




#### Make heatmap representations of CMat and PMat

# take away underscores
expabund$host.sp2 <- gsub("_"," " ,expabund$host.sp,)
expabund$par.sp2 <- gsub("_"," " ,expabund$par.sp,)

expabund$par.sp3 <- gsub("[a-z]*_",". ",expabund$par.sp)


# # order by summed abundance
expabundsort <- expabund %>% group_by(par.sp) %>% mutate(parsum = sum(beta))
expabundsort2 <- expabundsort %>% group_by(host.sp) %>% mutate(hostsum = sum(beta))

expabundsort2$par.sp <- gsub("_"," " ,expabundsort2$par.sp)
expabundsort2$host.sp <- gsub("_"," " ,expabundsort2$host.sp)

expabundsort2 <- expabundsort2[order(expabundsort2$parsum,expabundsort2$hostsum,decreasing=TRUE),]



expabundsort2$host.sp <- with(expabundsort2,factor(host.sp,levels = (unique(as.character(host.sp)))))
expabundsort2$par.sp <- with(expabundsort2,factor(par.sp,levels = rev(unique(as.character(par.sp)))))



base_size <- 12
pdf("~/dat/bipartitemodelsBC/figs/abundheatmap-sort-sp.pdf",height=12)
ggplot(expabundsort2, aes(host.sp,par.sp)) + geom_tile(aes(fill = beta), colour = "white") + scale_fill_gradient(low = "white", high = "black") + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.8, angle = 270, hjust = 0, colour = "grey50"),axis.text.y=element_text(size=base_size*0.7))
dev.off()



#base_size <- 12
pdf("~/dat/bipartitemodelsBC/figs/abundheatmap-abbrev-tall2.pdf",height=12)
ggplot(expabund, aes(host.sp2,par.sp3)) + geom_tile(aes(fill = beta), colour = "white") + scale_fill_gradient(low = "white", high = "black") + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.8, angle = 270, hjust = 0, colour = "grey50"),axis.text.y=element_text(size=base_size*0.7))
dev.off()

PMatmean <- t(apply(output$BUGSoutput$sims.list$use, c(2, 3), mean))
dimnames(PMatmean) <- list(levels(long$par.sp),levels(long$host.sp))
prob.use2<- melt(PMatmean)
prob.use2<-prob.use2[,c(2,1,3)]
names(prob.use2)<-c( 'host.sp', 'par.sp', 'use')

prob.use2$host.sp2 <- gsub("_"," " ,prob.use2$host.sp)
prob.use2$par.sp2 <- gsub("_"," " ,prob.use2$par.sp)

prob.use2$par.sp3 <- gsub("[a-z]*_",". ",prob.use2$par.sp)


# #ALTERNATIVE WITH MEDIAN
# PMatmedian <- t(apply(output$BUGSoutput$sims.list$use, c(2, 3), median))
# dimnames(PMatmedian) <- list(levels(long$par.sp),levels(long$host.sp))
# prob.use3<- melt(PMatmedian)
# prob.use3<-prob.use3[,c(2,1,3)]
# names(prob.use3)<-c( 'host.sp', 'par.sp', 'use')


pdf("~/dat/bipartitemodelsBC/figs/useheatmap-tall.pdf",height=12)
ggplot(prob.use2, aes(host.sp2,par.sp3)) + geom_tile(aes(fill = use), colour = "white") + scale_fill_gradient(low = "white", high = "black") + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.8, angle = 270, hjust = 0, colour = "grey50"),axis.text.y=element_text(size=base_size*0.7,vjust=0.25))
dev.off()


# # order by summed use
probusesort <- prob.use2 %>% group_by(par.sp) %>% mutate(parsum = sum(use))
probusesort2 <- probusesort %>% group_by(host.sp) %>% mutate(hostsum = sum(use))

probusesort2$par.sp <- gsub("_"," " ,probusesort2$par.sp)
probusesort2$host.sp <- gsub("_"," " ,probusesort2$host.sp)

probusesort2 <- probusesort2[order(probusesort2$parsum,probusesort2$hostsum,decreasing=FALSE),]

probusesort2$host.sp <- with(probusesort2,factor(host.sp,levels = rev(unique(as.character(host.sp)))))
probusesort2$par.sp <- with(probusesort2,factor(par.sp,levels = (unique(as.character(par.sp)))))

pdf("~/dat/bipartitemodelsBC/figs/useheatmap-sort2.pdf",height=12)
ggplot(probusesort2, aes(host.sp,par.sp)) + geom_tile(aes(fill = use), colour = "white") + scale_fill_gradient(low = "white", high = "black") + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.8, angle = 270, hjust = 0, colour = "grey50"),axis.text.y=element_text(size=base_size*0.7,vjust=0.25))

dev.off()

prob.use$presence <- ifelse(prob.use$prev > 0, 1,0)
hostuse <- prob.use %>% group_by(host.sp) %>% summarize(degree=sum(presence))
# https://johnbaumgartner.wordpress.com/2012/06/07/r-functions-to-filter-rjags-results/
library(jagstools)
ID<-as.data.frame(jagsresults(output, params='PD_host'))
names(ID)<-c("mean","sd","q2.5","q25","q50","q75","q97.5","Rhat","n.eff")
row.names(ID)<-levels(long$host.sp)
ID$degree <- hostuse$degree
ID<-ID[order(ID[,1]),]
ID$species <- gsub("_"," " ,row.names(ID))

pdf("~/dat/bipartitemodelsBC/figs/host-pardiv.pdf")
par(oma=c(4, 11, 0, 1), mar=c(0, 0.3, 2, 0), bg ='white') # try changing all margins to 0??
plot.new()
with(ID[1:nrow(ID),], {
  plot.window(xlim = range(c(ID$q2.5,ID$q97.5)), ylim =range(c(1, nrow(ID))))
  box()
  axis(1)
  axis(2, at=1:nrow(ID), las=2, labels=ID$species,cex.axis=1, tck =0.0)
  # title(main = "Insect diversity", line=1)
  title(xlab= "Parasite diversity", outer=TRUE)
  segments(x0=ID$q2.5, y0=c(1:nrow(ID)), x1=ID$q97.5, y1=c(1:nrow(ID)), lty=1, lwd=2.5)
  segments(x0=ID$q25, y0=c(1:nrow(ID)), x1=ID$q75, y1=c(1:nrow(ID)), lty=1, lwd=3.5)
  points(x=ID$q50, y=c(1:nrow(ID)), cex=1.5, lwd=1.5, pch=21, bg="red" )
  #abline(v=0, col=8)
  points(x=ID$degree,y=c(1:nrow(ID)),pch=4)
})
dev.off()


paruse <- prob.use %>% group_by(par.sp) %>% summarize(degree=sum(presence))
HB<-as.data.frame(jagsresults(output, params='HB_invert'))
names(HB)<-c("mean","sd","q2.5","q25","q50","q75","q97.5","Rhat","n.eff")
row.names(HB)<-levels(long$par.sp)

HB$species <- row.names(HB)
HB <- merge(HB,paruse,by.x="species",by.y="par.sp",all=TRUE)
HB<-HB[order(HB[,"mean"]),]
HB$species2 <- gsub("_"," " ,HB$species)

pdf("~/dat/bipartitemodelsBC/figs/par-hostbr.pdf",height=11)
par(oma=c(4, 11, 0, 1), mar=c(0, 0.3, 2, 0))
plot.new()
with(HB[1:nrow(HB),], {
  plot.window(xlim = range(c(HB$q2.5,HB$q97.5)), ylim =range(c(1, nrow(HB))))
  box()
  axis(1, at=c(0:10), labels= c(0:10), )
  axis(2, at=1:nrow(HB), labels=HB$species2, las=1, cex.axis=0.7, tck =0.0)
  #axis(2, at = 20, labels="Insect species", las=3, cex.axis=1, tck =0.0)
  #title(main = "Host Breadth", line=1)
  title(xlab= "Number of hosts", outer=TRUE)
  segments(x0=HB$q2.5, y0=c(1:nrow(HB)), x1=HB$q97.5, y1=c(1:nrow(HB)), lty=1, lwd=2.5)
  segments(x0=HB$q25, y0=c(1:nrow(HB)), x1=HB$q75, y1=c(1:nrow(HB)), lty=1, lwd=3.5)
  points(x=HB$q50, y=c(1:nrow(HB)), cex=1.2, lwd=1, pch=21, bg="red" )
  points(x=HB$degree,y=(1:nrow(HB)),pch=4)
})
dev.off()

agg$presence <- as.numeric(agg$npos > 0)

########################################


library(bipartite)
#  make matrices - higher level (par) in columns, lower level (host) in rows
colnames <- c("mean","sd","q2.5","q25","q50","q75","q97.5","Rhat","n.eff")

use.out <- as.data.frame(jagsresults(output,params="use"))
names(use.out) <- colnames
use.out$par.sp <- rep(levels(long$par.sp),each=long$Nhost.sp)
use.out$host.sp <- rep(levels(long$host.sp),long$Npar)

usemerge <- merge(agg,use.out,by=c("host.sp","par.sp"))

mat.obs <- usemerge %>% select(host.sp,par.sp,presence) %>% spread(par.sp,presence) %>% select(-host.sp) %>% as.matrix()
rownames(mat.obs) <- levels(usemerge$host.sp)

mat.low <- usemerge %>% select(host.sp,par.sp,q2.5) %>% spread(par.sp,q2.5) %>% select(-host.sp) %>% as.matrix()
rownames(mat.low) <- levels(usemerge$host.sp)

mat.med <- usemerge %>% select(host.sp,par.sp,q50) %>% spread(par.sp,q50) %>% select(-host.sp) %>% as.matrix()
rownames(mat.med) <- levels(usemerge$host.sp)

mat.hi <- usemerge %>% select(host.sp,par.sp,q97.5) %>% spread(par.sp,q97.5) %>% select(-host.sp) %>% as.matrix()
rownames(mat.hi) <- levels(usemerge$host.sp)

m.0 <- networklevel(mat.obs,index="binary")

m.l <- networklevel(mat.low,index="binary")
m.m <- networklevel(mat.med,index="binary") # topology
m.h <- networklevel(mat.hi,index="binary") # topology

all(mat.low==mat.obs)
all(mat.med==mat.obs)
sum(mat.hi!=mat.obs)

bip.all <- data.frame(m0=m.0,ml=m.l,mm=m.m,mh=m.h)

cbind(m.m,m.h)


library(jagstools)

host<-as.data.frame(jagsresults(output, params='hosts'))
names(host) <- colnames
#dim(output$BUGSoutput$sims.list$hosts) # this gives 3000*16*16 array
# extract row and col to convert to matrix
host$names <- rownames(host)
host$row <- substring(host$names,regexpr("\\[",host$names)+1,regexpr(",",host$names)-1)
host$col <- substring(host$names,regexpr(",",host$names)+1,regexpr("\\]",host$names)-1)

# observed
obshost <- mat.obs %*% t(mat.obs)

lowhost<- matrix(0,nrow=16,ncol=16,dimnames=dimnames(obshost))
for (i in 1:16) {
  for (j in 1:16) {
    ind <- which(host$row==i & host$col ==j)
    lowhost[i,j] <- host$q2.5[ind] 
  }
}


medianhost <- matrix(0,nrow=16,ncol=16,dimnames=dimnames(obshost))
for (i in 1:16) {
  for (j in 1:16) {
    ind <- which(host$row==i & host$col ==j)
    medianhost[i,j] <- host$q50[ind] 
  }
}

highhost <- matrix(0,nrow=16,ncol=16,dimnames=dimnames(obshost))
for (i in 1:16) {
  for (j in 1:16) {
    ind <- which(host$row==i & host$col ==j)
    highhost[i,j] <- host$q97.5[ind] 
  }
}

all(lowhost==medianhost)
which(lowhost!=obshost)

library(igraph)
g0 <- graph.adjacency(obshost,mode="undirected",diag=FALSE,weighted=TRUE)
g2 <- graph.adjacency(medianhost,mode="undirected",diag=FALSE,weighted=TRUE)
g1 <- graph.adjacency(lowhost,mode="undirected",diag=FALSE,weighted=TRUE)
g3 <- graph.adjacency(highhost,mode="undirected",diag=FALSE,weighted=TRUE)
# let's remove multi-edges and loops
#g <- simplify(g)
# g <- g1

# let's see if we have communities here using the 
# Grivan-Newman algorithm
# 1st we calculate the edge betweenness, merges, etc...
ebc0 <- edge.betweenness.community(g0, directed=F)
ebc1 <- edge.betweenness.community(g1, directed=F)
ebc2 <- edge.betweenness.community(g2, directed=F)
ebc3 <- edge.betweenness.community(g3, directed=F)

communities(ebc0)
communities(ebc1)
communities(ebc2)
communities(ebc3)

# Now we have the merges/splits and we need to calculate the modularity
# for each merge for this we'll use a function that for each edge
# removed will create a second graph, check for its membership and use
# that membership to calculate the modularity
# mods <- sapply(0:ecount(g0), function(i){
#   gTEMP <- delete.edges(g0, ebc0$removed.edges[seq(length=i)])
#   cl <- clusters(gTEMP)$membership
# # March 13, 2014 - compute modularity on the original graph g 
# # (Thank you to Augustin Luna for detecting this typo) and not on the induced one g2. 
#   modularity(g0,cl)
# })
#  
# # we can now plot all modularities
# plot(mods, pch=20)

# Now, let's color the nodes according to their membership
# gTEMP<-delete.edges(g0, ebc0$removed.edges[seq(length=which.max(mods)-1)])
# V(g0)$color=clusters(gTEMP)$membership
#  
# # Let's choose a layout for the graph
# g0$layout <- layout.fruchterman.reingold
#  
# # plot it
# plot(g2, vertex.label=NA)

# define colors
V(g0)$color <- membership(ebc0)
V(g1)$color <- membership(ebc1)
V(g2)$color <- membership(ebc2)
V(g3)$color <- membership(ebc3)

# if we wanted to use the fastgreedy.community agorithm we would do
fc0 <- fastgreedy.community(g0)
fc1 <- fastgreedy.community(g1)
fc2 <- fastgreedy.community(g2)
fc3 <- fastgreedy.community(g3)

communities(fc0)
fc0
communities(fc3)
fc3
communities(fc1)
fc1
communities(fc2)
fc2
# com<-community.to.membership(g0, fc$merges, steps= which.max(fc$modularity)-1)
# V(g)$color <- com$membership+1
# g$layout <- layout.fruchterman.reingold
# plot(g, vertex.label=NA)

pdf("~/dat/bipartitemodelsBC/figs/networks2.pdf",height=10,width=10)
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(0,0,2,0))
# plot.igraph(g0,layout=layout.fruchterman.reingold, edge.color="grey",edge.width=E(g0)$weight,main="Observed")
# plot.igraph(g1,layout=layout.fruchterman.reingold, edge.color="grey",edge.width=E(g1)$weight,main="Lower bound")
# plot.igraph(g2,layout=layout.fruchterman.reingold, edge.color="grey",edge.width=E(g2)$weight,main="Median")
# plot.igraph(g3,layout=layout.fruchterman.reingold, edge.color="grey",edge.width=E(g3)$weight,main="Upper bound")
# dev.off()

plot.igraph(g0,layout=layout.circle(g0), edge.color="grey",edge.width=E(g0)$weight,main="Observed")
plot.igraph(g1,layout=layout.circle(g1), edge.color="grey",edge.width=E(g1)$weight,main="Lower bound")
plot.igraph(g2,layout=layout.circle(g2), edge.color="grey",edge.width=E(g2)$weight,main="Median")
plot.igraph(g3,layout=layout.circle(g3), edge.color="grey",edge.width=E(g3)$weight,main="Upper bound")
dev.off()

# define colors
V(g0)$color <- membership(fc0)
V(g1)$color <- membership(fc1)
V(g2)$color <- membership(fc2)
V(g3)$color <- membership(fc3)


pdf("~/dat/bipartitemodelsBC/figs/networks3.pdf",height=10,width=10)
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(0,0,2,0))
plot.igraph(g0,layout=layout.circle(g0), edge.color="grey",edge.width=E(g0)$weight,main="Observed")
plot.igraph(g1,layout=layout.circle(g1), edge.color="grey",edge.width=E(g1)$weight,main="Lower bound")
plot.igraph(g2,layout=layout.circle(g2), edge.color="grey",edge.width=E(g2)$weight,main="Median")
plot.igraph(g3,layout=layout.circle(g3), edge.color="grey",edge.width=E(g3)$weight,main="Upper bound")
dev.off()


