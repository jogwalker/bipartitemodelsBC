# try to make a mammal phylogeny
#setwd("~/Documents/ch3 stuff")
# update 3 May 2016
# revise 7 July for thesis (using iTOL to make phylogeny)

library(ape)

tree <- read.tree("~/dat/bipartitemodelsBC/TimetreeOfLife2015.nwk")

species <- c("Connochaetes_taurinus",
             "Tragelaphus_scriptus",
             "Syncerus_caffer",
             "Sylvicapra_grimmia",
             "Oryx_gazella",
             "Giraffa_camelopardalis",
             "Tragelaphus_strepsiceros",
             "Aepyceros_melampus",
             "Alcelaphus_buselaphus",
             "Antidorcas_marsupialis",
             "Raphicerus_campestris",
             "Equus_burchellii",
             "Bos_taurus",
             "Bos_indicus",
             "Equus_asinus",
             "Equus_caballus",
             "Ovis_aries",
             "Hippopotamus_amphibius",
             "Ceratotherium_simum",
             "Loxodonta_africana",
             "Taurotragus_derbianus",
             "Hippotragus_equinus",
             "Hippotragus_niger",
             "Capra_hircus")

pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, species))
pruned.tree$tip.label <- c("Horse","Burchell's Zebra","Donkey","White Rhinoceros","Hippopotamus","Giraffe","Cattle (B. taurus)","Cattle (B. indicus)","Cape Buffalo","Greater Kudu","Bushbuck","Eland","Sable","Roan","Gemsbok","Red Hartebeest","Blue Wildebeest","Goat","Sheep","Springbok","Steenbok","Common Duiker","Impala","African Elephant")
write.tree(pruned.tree,file="~/git/bipartitemodelsBC/commontree.txt")
#http://itol.embl.de/upload.cgi



##################
########
times <- read.csv("~/dat/bipartitemodelsBC/time distance.csv",head=T)
library(tidyr)
library(dplyr)

timesl <- times %>% gather("Species2","time",2:ncol(times),na.rm=T)

library(rjags)
library(jagstools)
#load("~/Documents/bipartitemodelsBC/results/cnj_output-trunc1e+05.RData")
#load("~/Documents/bipartitemodelsBC/results/cnjlong-trunc1e+05.RData")
# already loaded march 300k from running finalplotcode 
host<-as.data.frame(jagsresults(output, params='hosts'))
host$names <- rownames(host)
host$row <- substring(host$names,regexpr("\\[",host$names)+1,regexpr(",",host$names)-1)
host$col <- substring(host$names,regexpr(",",host$names)+1,regexpr("\\]",host$names)-1)

# make species names comparable & merge
species <- levels(long$host.sp)
scispecies <- levels(timesl$Species)
sciorder <- c(5,15,14,4,13,6,10,9,16,7,1,2,11,3,12,8)
species <- data.frame(common=species,sci=scispecies[sciorder])
species$level <- 1:16
species$dot <- gsub(" ",".",species$sci)

timesl2 <- merge(timesl,species[2:3],by.x="Species",by.y="sci")
timesl3 <- merge(timesl2,species[3:4],by.x="Species2",by.y="dot")
names(timesl3)[4:5] <- c("col","row")

full <- merge(host,timesl3,by=c("col","row"))
names(full)[5:9] <- c("q2.5","q25","q50","q75","q97.5")

# add obshost
obshost <- data.frame(obshost)
obshost$species1.common <- row.names(obshost)
obslong <- gather(obshost,"species2.common","obscount",1:16)
obslong2 <- merge(obslong,species,by.x="species1.common",by.y="common")
obslong3 <- merge(obslong2,species,by.x="species2.common",by.y="common")
full2 <- left_join(full,obslong3,by=c("Species2"="dot.x","Species"="sci.y"))


plot(full$mean ~ full$time)

summary(lm(mean ~ time,data=full))
summary(lm(q50 ~ time,data=full))

cor.test(full2$mean,full2$time,method="spearman")
cor.test(full2$q97.5,full2$time,method="spearman")
cor.test(full2$obscount,full2$time,method="spearman")
cor.test(full2$q50,full2$time,method="spearman")
cor.test(full2$q2.5,full2$time,method="spearman")
cor.test(full2$q25,full2$time,method="spearman")
