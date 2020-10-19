setwd()
#loading libraries
library(rioja)
library(vegan)
library(ggplot2)
library(analogue)
library(velociraptr)
library(tidyverse)
library(grid)


#getting the data
pollen<-read.table("pollen.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
dictionary<-read.table("dictionary.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)


#defining time, depth and spike objects
depth<-pollen [, "depth"]
age<-pollen [, "kaBP"]
accrate<-pollen[,"accrate"]
lyc<-pollen [,"Lycopodium"]
lyc.conc<-pollen [,"Lyc.conc"]
undet <- pollen[,"Undetermined"]

# Estimating the sum of all pollen types, including spores
##########################################################
all.pollen.types <- dictionary [,"taxa"]
all.pollen.sum <-pollen [, all.pollen.types]

#percentages for upland pollen types
################################################
pollen.types<-dictionary[(dictionary$functional!="algae"
                          & dictionary$functional!="fern"
                          & dictionary$functional!="hydrophyte" 
                          & dictionary$functional!="hygrophyte"),"taxa"]

pollen.sum<-pollen[, pollen.types]


pollen.non.aquatics<-pollen.sum/rowSums(pollen.sum) * 100


#percentages for aquatic pollen types
################################################
aquatic.pollen.types<-dictionary[(dictionary$functional!="algae"
                                & dictionary$functional!="fern"
                                & dictionary$functional!="tree"
                                & dictionary$functional!="shrub"
                                & dictionary$functional!="herb") ,"taxa"]

pollen.aquatics<-pollen[, aquatic.pollen.types]

pollen.aquatics<-pollen.aquatics/rowSums(pollen.sum) * 100


#percentages for fern spores
################################################

spores.fern<-dictionary[dictionary$functional=="fern","taxa"]

fern.spores<-pollen[, spores.fern]

fern.spores <- fern.spores/rowSums(pollen.sum)*100

fern.spores <-data.frame(fern.spores)

#percentages for algae spores
################################################

spores.algae<-dictionary[dictionary$functional=="algae","taxa"]

algae.spores<-pollen[, spores.algae]

algae.spores <- algae.spores/rowSums(all.pollen.sum)*100

algae.spores <-data.frame(algae.spores)

##Percentages for functional groups
################################################

#getting the functional groups
functional.groups<-unique(dictionary$functional)

#removing algae
functional.groups<-functional.groups[functional.groups!="algae"]

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

#removing taxa in non bio-group
bio.groups<-bio.groups[bio.groups!="Botryococcus"] 
bio.groups<-bio.groups[bio.groups!="Cerealia"]
bio.groups<-bio.groups[bio.groups!="Chorcorus"]
bio.groups<-bio.groups[bio.groups!="melastomataceae"]
bio.groups<-bio.groups[bio.groups!="Poaceae"]
bio.groups<-bio.groups[bio.groups!="Aster"]
bio.groups<-bio.groups[bio.groups!="ruderal"]

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

save(bio.groups.df, file = "rawBioGroups.RData")

#computing percentage
bio.groups.df<-bio.groups.df/rowSums(bio.groups.df) * 100


##Binding all tables together
#############################################3

pollen.p<-cbind(bio.groups.df,functional.groups.df, pollen.non.aquatics, pollen.aquatics, algae.spores, fern.spores)
pollen.p<-cbind(age, depth, pollen.p)

#writing result
save(pollen, file = "pollen.RData")
write.table(pollen.p, file = "pollen_p.csv", row.names = FALSE, col.names = TRUE, sep=",")


# Calculating PAR and influxes
#####################################################
#####################################################

#create all objects needed to estimate individual PAR
pollen.types<-dictionary[dictionary$functional!="algae","taxa"] #selecting taxa to PAR

pollen.par<-pollen[, pollen.types]#assigning values to them

total.pollen <-rowSums(pollen.par)#total pollen sum 

names(pollen.par)

#create table to populate with the PAR results
pollen.par.df<-data.frame(matrix(NA, ncol=ncol(pollen.par), nrow=nrow(pollen.par)))
colnames(pollen.par.df)<-names(pollen.par)

#loop and populate dataframe
for (pollen.par.df in pollen.par.df){
  
  pollen.par.df<-((pollen.par/lyc)*lyc.conc)/accrate
  
}

#calculate PAR for bioclimatic groups 
load("rawBioGroups.RData")
#create table to populate with the PAR bioclimatic groups results
bio.par.df<-data.frame(matrix(NA, ncol=ncol(bio.groups.df), nrow=nrow(bio.groups.df)))
colnames(bio.par.df)<-names(bio.groups.df)
#loop and populate dataframe for bioclim groups
for (bio.par.df in bio.par.df){
  
  bio.par.df<-((bio.groups.df/lyc)*lyc.conc)/accrate
  
}

#calculate PAR for the total pollen sum
total.pollen.par <-((total.pollen/lyc)* lyc.conc)/accrate

#merge total pollen PAR and taxa PAR
pollen.par <-cbind(pollen.par.df,total.pollen.par)

write.table(pollen.par, file="pollen_par.csv", row.names = FALSE, col.names = TRUE, sep=",")

#save bioclim groups PAR for plotting later
save(bio.par.df, file = "bioPAR.RData")


##Estimating principal curves for pollen percentages
################################################################################
################################################################################
pollen_p <- read.csv("pollen_p.csv", header=TRUE, stringsAsFactors = FALSE)

pollen.p <- pollen_p[-1:-2]#remove age and depth

#In order to test different Principal curves that might be indicating 
#individual or collective gradients, I chose to create PrC for various 
#taxa groups.
#######################################################################

prcurve.bio <- pollen_p[, c("Lower.forest.limit", "Forest.belt", "Afroalpine",
                            "Ericaceous.belt","Upper.forest.limit")]

# Change any NA it may have changed from 0 back to 0
prcurve.bio[is.na(prcurve.bio)] = 0

# Remove taxa where total abundanace is less than 2% try both 
mx <- apply(prcurve.bio, 2, max)
prcurve.bio2 <- prcurve.bio[, mx>2]

#Principal Curve curve fitting through CA and PCA
#CA increasing the number of iterations to 50
p.pcca <- prcurve(prcurve.bio2, method = "ca", maxit=50, 
                  trace = T, vary = T, penalty = 1.4, thresh = 0.0005, plotit = TRUE)#CA 
p.pcca
summary(p.pcca)
plot(p.pcca)

#PCA 
p.pcpc <- prcurve(prcurve.bio2, method = "pca", maxit=50, 
                  thres= 5e-04, trace = TRUE, vary = FALSE, penalty = 1.4, plotit=T)#PCA
p.pcpc
summary(p.pcpc)
plot(p.pcpc)

#Compare explained variance with axes of other tests 
pca <- varExpl(rda(prcurve.bio2), axes=1:4, cumulative = T)
ca <- varExpl(cca(prcurve.bio2), axes=1:4, cumulative = T)
pca
ca

## Extract position of the Principal Curve scores in the produced object
pos.bio <- scores(p.pcpc, display = "curve")
head(pos.bio)


##Saving the Principal Curve for bioclimatic groups scores for plotting
write.table(pos.bio, file="PrC_bio.csv", row.names=FALSE, col.names=TRUE, sep=",")

##Creating a PrC for individual taxa
#######################################
prcurve.p <- pollen_p[-1:-13]#remove functional groups and aquatics
prcurve.p <- prcurve.p[-91:-98]
# Change any NA it may have changed from 0 back to 0
prcurve.p[is.na(prcurve.p)] = 0
# Remove taxa where total abundanace is less than 2% try both 
mx <- apply(prcurve.p, 2, max)
prcurve.p2 <- prcurve.p[, mx>2]

#Principal Curve curve fitting through CA and PCA
#CA increasing the number of iterations to 50
p.pcca <- prcurve(prcurve.p2, method = "ca", maxit=50, 
                  trace = T, vary = T, penalty = 1.4, thresh = 0.0005, plotit = TRUE)#CA 
p.pcca
summary(p.pcca)
plot(p.pcca)

#PCA 
p.pcpc <- prcurve(prcurve.p2, method = "pca", maxit=50, 
                  thres= 5e-04, trace = TRUE, vary = FALSE, penalty = 1.4, plotit=T)#PCA
p.pcpc
summary(p.pcpc)
plot(p.pcpc)

#Compare explained variance with axes of other tests 
pca <- varExpl(rda(prcurve.p2), axes=1:4, cumulative = T)
ca <- varExpl(cca(prcurve.p2), axes=1:4, cumulative = T)
pca
ca
## Extract position of the Principal Curve scores in the produced object
pos.p2 <- scores(p.pcpc, display = "curve")
head(pos.p2)

##Saving the Principal Curve scores for plotting
write.table(pos.p2, file="PrC_p2.csv", row.names=FALSE, col.names=TRUE, sep=",")

#Joining Principal Curve scores with pollen data, creating RData
################################################################

pollen.p<-cbind(pollen.p, pos.p2, pos.bio)
save(pollen.p, file = "pollen.p.RData")

##Plotting percentage pollen diagram
#################################################################################
#################################################################################

# Select taxa above the 2%
mx <- apply(pollen.p, 2, max)
pollen.p2 <- pollen.p[, mx>2]

## Change any NA it may have changed from 0 back to 0
pollen.p2[is.na(pollen.p2)] = 0

#create table with pecentages >2%
write.table(pollen.p2, file="pollen.p2.csv", sep= ",")

##plotting diagram using strat.plot
ylim <- range (0,13.700)
ymin <-0.094
ymax <-13.700
y.scale.1 <-seq (0, 0.094)
y.scale.2<-seq (1, 14, by=0.5)
y.tks<- sort(c(y.scale.1,y.scale.2))

#total diagram
names(pollen.p2)
dev.off()
pollen.diagram <- strat.plot(d=pollen.p2, scale.percent=TRUE, yvar=age, y.rev=TRUE,
                          ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.bar=FALSE, 
                          plot.line = FALSE,x.pc.lab=TRUE, col.poly="forestgreen")

#Clustering and pollen zones
#####################################################
dissc <- dist(sqrt(pollen.sum/100)^2)
clust <- chclust(dissc, method = "coniss")
# broken stick model suggests significant zones and add to diagram
bstick(clust)
plot(clust)#plotting the dendrogram
plot(clust, hang=-1, horiz=TRUE)# Rotated through 90 degrees

# Rotated and observations plotted according to sample depth.
plot(clust, xvar=pollen_p$age, hang=-1, horiz=TRUE, x.rev=TRUE)
addClustZone(pollen.diagram, clust,6, col="red")

##Plotting PAR diagram
#################################################################################
#################################################################################
load("bioPAR.RData")

dev.off()
bioPAR.diagram <- strat.plot(d=bio.par.df, scale.percent=FALSE, x.pc.omit0=FALSE, yvar=age, y.rev=TRUE,
                             ylim=ylim, y.tks=y.tks, plot.poly=FALSE, plot.bar=TRUE,
                             lwd.bar=2, plot.line = FALSE, col.bar="forestgreen")

