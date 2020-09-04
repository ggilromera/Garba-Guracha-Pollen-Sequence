setwd()
#loading libraries
library(rioja)
#library(ggpalaeo)
library(vegan)
library(ggplot2)
library(analogue)
library(velociraptr)
library(tidyverse)
library(grid)



#getting the data
pollen<-read.table("pollen_v5.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
dictionary<-read.table("dictionary_v6.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

#percentages for upland pollen types
################################################
pollen.types<-dictionary[(dictionary$functional!="algae"
                          & dictionary$functional!="fern"
                          & dictionary$functional!="hydrophyte" 
                          & dictionary$functional!="hygrophyte"),"taxa"]

pollen.non.aquatics<-pollen[, pollen.types]

pollen.non.aquatics<-pollen.non.aquatics/rowSums(pollen.non.aquatics) * 100


#percentages for aquatic pollen types
################################################
aquatic.pollen.types<-dictionary[(dictionary$functional=="hydrophyte"
                                & dictionary$functional=="hygrophyte") ,"taxa"]

pollen.aquatics<-pollen[, aquatic.pollen.types]

#this is the total pollen assemblage on which the percentages of non tree, shrubs, herbs are based
pollen.total <- dictionary[(dictionary$functional!="algae" 
                            & dictionary$functional!="fern"),"taxa"]

pollen.total <- pollen [, pollen.total] 


pollen.aquatics<-pollen.aquatics/rowSums(pollen.total) * 100

#percentages for algae spores
################################################

spores.algae<-dictionary[dictionary$functional=="algae","taxa"]

algae.spores<-pollen[, spores.algae]

algae.spores<-algae.spores/rowSums(pollen.total) * 100

algae.spores <-data.frame(algae.spores)

#percentages for faecal spores
################################################

spores.faecal<-dictionary[dictionary$functional=="fern","taxa"]

faecal.spores<-pollen[, spores.faecal]

pollen.total.fern <- dictionary[dictionary$functional!="algae","taxa"]

pollen.total.fungi <- pollen [, pollen.total.fungi]

faecal.spores<-faecal.spores/rowSums(pollen.total.fungi) * 100


##Percentages for functional groups
################################################

#getting the functional groups
functional.groups<-unique(dictionary$functional)

#removing algae
functional.groups<-functional.groups[functional.groups!="algae"]

#removing aquatics
functional.groups<-functional.groups[functional.groups!="aquatic.herb"]

#removing faecal spores
functional.groups<-functional.groups[functional.groups!="fungi"]

#create empty dataframe for the functional groups
functional.groups.df<-data.frame(matrix(NA, ncol=length(functional.groups), nrow=nrow(pollen)))
names(functional.groups.df)<-functional.groups

for (functional.group in functional.groups){
  
  #selecting the columns in pollen belonging to functional.group
  pollen.functional.group<-pollen[, dictionary[dictionary$functional==functional.group ,"taxa"]]
  
  #summing the rows (if more than one column, we use rowSums, otherwise we get the column itself)
  if(ncol(pollen.functional.group)>1){
    functional.groups.df[,functional.group]<-rowSums(pollen.functional.group)
  } else {
    functional.groups.df[,functional.group]<-pollen.functional.group
    }
  
}

#computing percentage
functional.groups.df<-functional.groups.df/rowSums(functional.groups.df) * 100

##Percentages for bioclimatic groups
################################################
#getting the bio groups
bio.groups<-unique(dictionary$bio)

#removing algae
bio.groups<-bio.groups[bio.groups!="Botryococcus"]

#removing aquatics
bio.groups<-bio.groups[bio.groups!="Cerealia"]

#removing faecal spores
bio.groups<-bio.groups[bio.groups!="Chorcorus"]
melastomataceae
Poaceae
"

#create empty dataframe for the bio groups
bio.groups.df<-data.frame(matrix(NA, ncol=length(bio.groups), nrow=nrow(pollen)))
names(bio.groups.df)<-bio.groups

for (bio.group in bio.groups){
  
  #selecting the columns in pollen belonging to bio.group
  pollen.bio.group<-pollen[, dictionary[dictionary$bio==bio.group ,"taxa"]]
  
  #summing the rows (if more than one column, we use rowSums, otherwise we get the column itself)
  if(ncol(pollen.bio.group)>1){
    bio.groups.df[,bio.group]<-rowSums(pollen.bio.group)
  } else {
    bio.groups.df[,bio.group]<-pollen.bio.group
  }
  
}

#computing percentage
bio.groups.df<-bio.groups.df/rowSums(bio.groups.df) * 100


##Binding all tables together
#############################################3

pollen<-cbind(pollen, functional.groups.df, pollen.non.aquatics, pollen.aquatics, algae.spores, faecal.spores)

age <-pollen$kaBP
depth <-pollen$depth

pollen.p<-cbind(age, depth, functional.groups.df, pollen.non.aquatics, pollen.aquatics, algae.spores, faecal.spores)

#writing result
save(pollen, file = "pollen.RData")
write.table(pollen.p, file = "pollen_p.csv", row.names = FALSE, col.names = TRUE, sep=",")


##Estimating principal curves for pollen percentages
################################################################################
################################################################################
pollen_p <- read.csv("pollen_p.csv", header=TRUE, stringsAsFactors = FALSE)

pollen.p <- pollen_p[-1:-2]#remove age and depth
names(pollen.p)

# Change any NA it may have changed from 0 back to 0
pollen.p[is.na(pollen.p)] = 0

# Remove taxa where total abundanace is less than 2% try both 
mx <- apply(pollen.p, 2, max)
pollen.p2 <- pollen.p[, mx>2]

#Principal Curve curve fitting through CA and PCA
#CA increasing the number of iterations to 50
p.pcca <- prcurve(pollen.p2, method = "ca", maxit=50, 
                  trace = T, vary = T, penalty = 1.4, thresh = 0.0005, plotit = TRUE)#CA 
p.pcca
summary(p.pcca)
plot(p.pcca)

#PCA 
p.pcpc <- prcurve(pollen.p2, method = "pca", maxit=50, 
                  thres= 5e-04, trace = TRUE, vary = FALSE, penalty = 1.4, plotit=T)#PCA
p.pcpc
summary(p.pcpc)
plot(p.pcpc)

## Extract position of the Principal Curve scores in the produced object
pos <- scores(p.pcca, display = "curve")
head(pos)

#Compare explained variance with axes of other tests 
pca <- varExpl(rda(pollen.p2), axes=1:4, cumulative = T)
ca <- varExpl(cca(pollen.p2), axes=1:4, cumulative = T)
pca
ca


##Saving the Principal Curve scores for plotting
write.table(pos, file="PrC.csv", row.names=FALSE, col.names=TRUE, sep=",")

#Joining Principal Curve scores with pollen data 
pollen_p_pc<-cbind(pollen_p, pos)

##Saving the Principal Curve scores for plotting
write.table(pollen_p_pc, file="pollen_p_pc.csv", row.names=FALSE, col.names=TRUE, sep=",")


##Plotting pollen diagram
#################################################################################
#################################################################################

## Change any NA it may have changed from 0 back to 0
pollen.p2[is.na(pollen.p2)] = 0

#create table with pecentages >2%
write.table(pollen.p2, file="pollen.p2.csv", sep= ",")

##plotting diagram using strat.plot
ylim <- range (-0.066,13.700)
ymin <--0.066
ymax <-13.700
y.scale.1 <-seq (-0.066,0, by=0.066)
y.scale.2<-seq (1, 14, by=1)
y.tks<- sort(c(y.scale.1,y.scale.2))

#total diagram
names(pollen.p2)
dev.off()
pollen.diagram <- strat.plot(d=pollen.p2, scale.percent=TRUE, yvar=age, y.rev=TRUE, 
                             ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.bar=FALSE, plot.line = FALSE, col.poly="forestgreen", title="Garba Guracha (3950m asl")

#partial diagrams
names(pollen.p2)
all.trees <-pollen.p2 [5:14]
names(all.trees)
p.col.tree <- c(rep("forestgreen", times=1), rep("darkgoldenrod3", times=2), rep("darkolivegreen2", times=7))
trees <-strat.plot(d=all.trees, scale.percent=TRUE, yvar=age, y.rev=TRUE, 
                   ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.bar = FALSE, plot.line = FALSE, col.poly = p.col.tree, col.line = p.col.tree,
                   title="Garba Guracha - tree (3950m asl)")


all.shrubs.herbs <-pollen.p2 [15:32]
names(all.shrubs.herbs)
p.col.shrubs.herbs <- c(rep("darkgoldenrod3", times=7), rep("wheat4", times=1), rep("darkolivegreen2", times=10))
ex <- c(rep(TRUE, times=2), FALSE, FALSE, rep(TRUE, times=6), FALSE, FALSE, rep(TRUE, times=6))
shrubs.herbs <- strat.plot(d=all.shrubs.herbs, scale.percent = TRUE, yvar=age, y.rev=TRUE, 
                           ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.bar=FALSE, plot.line = FALSE, col.poly=p.col.shrubs.herbs, 
                           exag=ex, col.exag="auto", exag.alpha=0.50, exag.mult=3,title="Garba Guracha - shrubs and herbs (3950m asl)")

aquatics.faecal <- pollen.p2 [33:37]
names(aquatics.faecal)
p.col.aq.fae <- c(rep("turquoise3", times=2), rep("slategray2", times=4))
aquatics.fae <-strat.plot(d=aquatics.faecal, scale.percent = TRUE, yvar=age, y.rev=TRUE, 
                      ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.bar=FALSE, plot.line = FALSE, col.poly  = p.col.aq.fae,
                      title="Garba Guracha - aquatics and faecal (3950m asl)")



#Clustering and pollen zones
#####################################################
dissc <- dist(sqrt(pollen.total/100)^2)
clust <- chclust(dissc, method = "coniss")
# broken stick model suggests significant zones and add to diagram
bstick(clust)
plot(clust)#plotting the dendrogram
plot(clust, hang=-1, horiz=TRUE)# Rotated through 90 degrees

# Rotated and observations plotted according to sample depth.
plot(clust, xvar=pollen_p$age, hang=-1, horiz=TRUE, x.rev=TRUE)
addClustZone(pollen.diagram, clust,6, col="red")

# Calculating and plotting PAR
#####################################################
pollen_conc<- read.csv("pollen_v5.csv", header=TRUE, stringsAsFactors = FALSE)
dictionary<- read.csv("dictionary.csv", header=TRUE, stringsAsFactors = FALSE)

#create all objects needed to estimate PAR
pollen.types<-dictionary[dictionary$functional!="algae","taxa"] #selecting taxa to PAR

pollen.par<-pollen_conc[, pollen.types]#assigning values to them

total.pollen <-rowSums(pollen.par)#total pollen estimation (all taxa but algae)

accrate <-pollen_conc$accrate
lyc <-pollen_conc$Lycopodium
names(pollen.par)

#create table to populate with the PAR results
pollen.par.df<-data.frame(matrix(NA, ncol=ncol(pollen.par), nrow=nrow(pollen.par)))
colnames(pollen.par.df)<-names(pollen.par)

#loop and populate dataframe
for (pollen.par.df in pollen.par.df){
  
    pollen.par.df<-((pollen.par/lyc)* 31272)/accrate
  
}

#calculate PAR for the total pollen sum
total.pollen.par <-((total.pollen/lyc)* 31272)/accrate

#merge total pollen PAR and taxa PAR
pollen.par <-cbind(pollen.par.df,total.pollen.par)

write.table(pollen.par, file="pollen_par.csv", row.names = FALSE, col.names = TRUE, sep=",")


## Create age and depth vectors
age <-pollen_conc$kaBP
depth <-pollen_conc$depth
accrate <-pollen_conc$accrate
##create age axis limits and tick marks
ylim <- range (-0.066,13.700)
ymin <--0.066
ymax <-13.700
y.scale.1 <-seq (-0.066,0, by=0.066)
y.scale.2<-seq (1, 14, by=1)
yticks<- sort(c(y.scale.1,y.scale.2))

##plotting diagram using Stratiplot from (analogue)

#partial diagrams
names(pollen.par)
dev.off()
pollen.par.spores <-pollen.par[1:4]
pollen.par.trees <-pollen.par[5:10]

names(pollen.par.spores)
graph.widths <- c(15000, 15000, 15000, 15000)
pollen.diagram.spores <- strat.plot(d=pollen.par.spores, scale.percent = FALSE, graph.widths = graph.widths, x.pc.inc=15, yvar=age, y.rev=TRUE, 
                             ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.bar=FALSE, plot.line = FALSE, title="Garba Guracha (3950m asl")

Stratiplot(pollen.par.spores, 
           y=age,
           ylim=ylim, 
           yticks=yticks,
           varTypes="absolute",
           sort = "wa", 
           type = "poly",
           ylab ='ka BP')

Stratiplot(pollen.par.trees, 
           y=age, 
           ylim=ylim, 
           varTypes="absolute",  
           sort = "wa", 
           type = "poly",
           ylab ='Years Before Present')
