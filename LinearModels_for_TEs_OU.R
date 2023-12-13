library(ggplot2)
library(ggdendro)
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(tidyverse)
library(ggtree)
library(scico)


#File Paths
TEdata_Path<-"~/Documents/PolypterusBichir/TEproject/TEdataWithLifeHistory.csv"
fishtree_Path<-"~/Documents/PolypterusBichir/TEproject/NoNodeLabelsUpdatedTimeTree.phy"

#Read in data
TEtypesdata <- read.csv(TEdata_Path, header = TRUE, row.names=1)

#data curation
TEtypesdata$taxonomy <- "teleost"
TEtypesdata["Polypterus_bichir_lapradei",]$taxonomy <- "non-teleost"
TEtypesdata["Polypterus_senegalus",]$taxonomy <- "non-teleost"
TEtypesdata["Acipenser_ruthenus",]$taxonomy <- "non-teleost"
#don't forget the reedfish!
TEtypesdata["Erpetoichthys_calabaricus",]$taxonomy <- "non-teleost"
TEtypesdata["Lepisosteus_osseus",]$taxonomy <- "non-teleost"
TEtypesdata["Lepisosteus_oculatus",]$taxonomy <- "non-teleost"
TEtypesdata["Amia_calva",]$taxonomy <- "non-teleost"

##-- scatter plot with the TE content percent and the genome size data
plot(log(TEtypesdata$GenomeSize), TEtypesdata$Total_TE, col = factor(TEtypesdata$taxonomy), pch=19)
barplot(TEtypesdata$GenomeSize, col = factor(TEtypesdata$taxonomy))

##--read in the tree and modify the tip labels
fishtree <- read.tree(fishtree_Path)
fishtree$tip.label[fishtree$tip.label=="Monocentris_japonicus"] <-"Monocentris_japonica"
fishtree$tip.label[fishtree$tip.label=="Haplochromis_burtoni"] <-"Astatotilapia_burtoni"
fishtree$tip.label[fishtree$tip.label=="Maylandia_zebra"] <-"Metriaclima_zebra"
fishtree$tip.label[fishtree$tip.label=="Polypterus_bichir"] <-"Polypterus_bichir_lapradei"


## CHECK THAT THE NAMES ALL MATCH BETWEEN TREE AND DATA
treedata(fishtree, TEtypesdata) #this function is in the geiger package
#this matches 

##--scale the tree to avoid errors. fix it by making a temp tree
tempTree <- fishtree
tempTree$edge.length <- tempTree$edge.length * 100

#make brownian correlation structure (do this first, then use it in the gls function)
cmOU<-corMartins(1, tempTree, form= ~Species)


##--pgls with total TE
pglsModel <- gls(Total_TE ~ GenomeSize, correlation = cmOU,
                 data = TEtypesdata, method = "ML")
summary(pglsModel)
plot(Total_TE~GenomeSize, data= TEtypesdata, col = factor(TEtypesdata $taxonomy), pch=19)

##--pgls with DNA
pglsModelDNA <- gls(DNA ~ GenomeSize, correlation = cmOU,
                    data = TEtypesdata, method = "ML")
summary(pglsModelDNA)
coef(pglsModelDNA)

##--pgls with LINE
pglsModelLINE <- gls(LINE ~ GenomeSize, correlation = cmOU,
                     data = TEtypesdata, method = "ML")
summary(pglsModelLINE)
coef(pglsModelLINE)

##pgls with SINE
pglsModelSINE <- gls(SINE ~ GenomeSize, correlation = cmOU,
                     data = TEtypesdata, method = "ML")
summary(pglsModelSINE)
coef(pglsModelSINE)

##--pgls with LTR
pglsModelLTR <- gls(LTR ~ GenomeSize, correlation = cmOU,
                    data = TEtypesdata, method = "ML")
summary(pglsModelLTR)
coef(pglsModelLTR)

#make one plot in ggplot
data_plot<-TEtypesdata%>%gather(TeType, Percent, 1:6)%>%filter(TeType!="Unclassified")

ggplot(data_plot, aes(x = GenomeSize, y = Percent, color = TeType, shape= taxonomy)) +
  geom_point(alpha=0.6) +
  labs(x = "GenomeLength", y = "Percent", color = "TE Type", shape = "Taxonomy") +
  scale_color_manual(values = c("#c91367", "#932598", "#fce2fb", "#fef5ee", "#380a8a")) +
  scale_shape_manual(values = c(17, 16)) +
  theme_classic()+
  geom_abline(intercept=coef(pglsModelDNA)[1], slope= coef(pglsModelDNA)[2], linetype = "dashed", col="#c91367") +
  geom_abline(intercept=coef(pglsModelLINE)[1], slope= coef(pglsModelLINE)[2], linetype = "dashed", col="#932598") +
  geom_abline(intercept=coef(pglsModelLTR)[1], slope= coef(pglsModelLTR)[2], linetype = "dashed", col="#fce2fb") +
  geom_abline(intercept=coef(pglsModelSINE)[1], slope= coef(pglsModelSINE)[2], linetype = "dashed", col="#fef5ee") +
  geom_abline(intercept=coef(pglsModel)[1], slope= coef(pglsModel)[2], linetype = "dashed", col="#380a8a") 

###################################################
####TE %s versus absolute latitude 
##################################################

#sampling is 100% so we can skip the name check and use the same cmOU object
##--pgls with total TE
pglsModel <- gls(Total_TE ~ abs(medianLat), correlation = cmOU,
                 data = TEtypesdata, method = "ML")
summary(pglsModel)
plot(Total_TE~ medianLat, data= TEtypesdata, col = factor(TEtypesdata $taxonomy), pch=19)

##--pgls with DNA
pglsModelDNA <- gls(DNA ~ abs(medianLat), correlation = cmOU,
                    data = TEtypesdata, method = "ML")
summary(pglsModelDNA)
coef(pglsModelDNA)

##--pgls with LINE
pglsModelLINE <- gls(LINE ~ abs(medianLat), correlation = cmOU,
                     data = TEtypesdata, method = "ML")
summary(pglsModelLINE)
coef(pglsModelLINE)

##pgls with SINE
pglsModelSINE <- gls(SINE ~ abs(medianLat), correlation = cmOU,
                     data = TEtypesdata, method = "ML")
summary(pglsModelSINE)
coef(pglsModelSINE)
plot(SINE ~ abs(medianLat), data= TEtypesdata, col = factor(TEtypesdata $taxonomy), pch=19)
abline(a = coef(pglsModelSINE)["(Intercept)"], b = coef(pglsModelSINE)["abs(medianLat)"], col = "red", lwd = 2)

##--pgls with LTR
pglsModelLTR <- gls(LTR ~ abs(medianLat), correlation = cmOU,
                    data = TEtypesdata, method = "ML")
summary(pglsModelLTR)
coef(pglsModelLTR)

#make one plot in ggplot
data_plot<-TEtypesdata%>%gather(TeType, Percent, 1:6)%>%filter(TeType!="Unclassified")

ggplot(data_plot, aes(x = abs(medianLat), y = Percent, color = TeType, shape= taxonomy)) +
  geom_point(alpha=0.6) +
  labs(x = "Median Latitude", y = "Percent", color = "TE Type", shape = "Taxonomy") +
  scale_color_manual(values = c("#c91367", "#932598", "#fce2fb", "#fef5ee", "#380a8a")) +
  scale_shape_manual(values = c(17, 16)) +
  theme_classic()+
  geom_abline(intercept=coef(pglsModelDNA)[1], slope= coef(pglsModelDNA)[2], linetype = "dashed", col="#c91367") +
  geom_abline(intercept=coef(pglsModelLINE)[1], slope= coef(pglsModelLINE)[2], linetype = "dashed", col="#932598") +
  geom_abline(intercept=coef(pglsModelLTR)[1], slope= coef(pglsModelLTR)[2], linetype = "dashed", col="#fce2fb") +
  geom_abline(intercept=coef(pglsModelSINE)[1], slope= coef(pglsModelSINE)[2], linetype = "dashed", col="#fef5ee") +
  geom_abline(intercept=coef(pglsModel)[1], slope= coef(pglsModel)[2], linetype = "dashed", col="#380a8a") 

###################################################
####TE %s versus log(bodysize) 
##################################################
#sampling is 100% so we can skip the name check and use the same cmOU object
##--pgls with total TE
pglsModel <- gls(Total_TE ~ log(size), correlation = cmOU,
                 data = TEtypesdata, method = "ML")
summary(pglsModel)
plot(Total_TE~ log(size), data= TEtypesdata, col = factor(TEtypesdata $taxonomy), pch=19)

##--pgls with DNA
pglsModelDNA <- gls(DNA ~ log(size), correlation = cmOU,
                    data = TEtypesdata, method = "ML")
summary(pglsModelDNA)
coef(pglsModelDNA)

##--pgls with LINE
pglsModelLINE <- gls(LINE ~ log(size), correlation = cmOU,
                     data = TEtypesdata, method = "ML")
summary(pglsModelLINE)
coef(pglsModelLINE)

##pgls with SINE
pglsModelSINE <- gls(SINE ~ log(size), correlation = cmOU,
                     data = TEtypesdata, method = "ML")
summary(pglsModelSINE)
coef(pglsModelSINE)
plot(SINE ~ abs(medianLat), data= TEtypesdata, col = factor(TEtypesdata $taxonomy), pch=19)

##--pgls with LTR
pglsModelLTR <- gls(LTR ~ log(size), correlation = cmOU,
                    data = TEtypesdata, method = "ML")
summary(pglsModelLTR)
coef(pglsModelLTR)

#make one plot in ggplot
data_plot<-TEtypesdata%>%gather(TeType, Percent, 1:6)%>%filter(TeType!="Unclassified")

ggplot(data_plot, aes(x = log(size), y = Percent, color = TeType, shape= taxonomy)) +
  geom_point(alpha=0.6) +
  labs(x = "log Body Size", y = "Percent", color = "TE Type", shape = "Taxonomy") +
  scale_color_manual(values = c("#c91367", "#932598", "#fce2fb", "#fef5ee", "#380a8a")) +
  scale_shape_manual(values = c(17, 16)) +
  theme_classic()+
  geom_abline(intercept=coef(pglsModelDNA)[1], slope= coef(pglsModelDNA)[2], linetype = "dashed", col="#c91367") +
  geom_abline(intercept=coef(pglsModelLINE)[1], slope= coef(pglsModelLINE)[2], linetype = "dashed", col="#932598") +
  geom_abline(intercept=coef(pglsModelLTR)[1], slope= coef(pglsModelLTR)[2], linetype = "dashed", col="#fce2fb") +
  geom_abline(intercept=coef(pglsModelSINE)[1], slope= coef(pglsModelSINE)[2], linetype = "dashed", col="#fef5ee") +
  geom_abline(intercept=coef(pglsModel)[1], slope= coef(pglsModel)[2], linetype = "dashed", col="#380a8a") 


###################################################
####TE %s versus depth
##################################################
data2<-TEtypesdata %>% na.omit()

## CHECK THAT THE NAMES ALL MATCH BETWEEN TREE AND DATA
td<-treedata(tempTree, data2) 
#use this object below to account for no depth data in non marine species 

#convert to data frame and make numeric
data.depth<-as.data.frame(td$data) %>% mutate_at(c(1:6,8:9), ~as.numeric(.))

#make brownian correlation structure (do this first, then use it in the gls function)
cmOU<-corMartins(1, td$phy, form= ~Species)

pglsModel <- gls(Total_TE ~ depth, correlation = cmOU,
                 data = data.depth, method = "ML")
summary(pglsModel)
plot(Total_TE~ depth, data= data.depth, col = factor(data.depth$taxonomy), pch=19)

##--pgls with DNA
pglsModelDNA <- gls(DNA ~ depth, correlation = cmOU,
                    data = data.depth, method = "ML")
summary(pglsModelDNA)
coef(pglsModelDNA)

##--pgls with LINE
pglsModelLINE <- gls(LINE ~ depth, correlation = cmOU,
                     data = data.depth, method = "ML")
summary(pglsModelLINE)
coef(pglsModelLINE)

##pgls with SINE
pglsModelSINE <- gls(SINE ~ depth, correlation = cmOU,
                     data = data.depth, method = "ML")
summary(pglsModelSINE)
coef(pglsModelSINE)

##--pgls with LTR
pglsModelLTR <- gls(LTR ~ depth, correlation = cmOU,
                    data = data.depth, method = "ML")
summary(pglsModelLTR)
coef(pglsModelLTR)

#make one plot in ggplot
data_plot<-data.depth%>%gather(TeType, Percent, 1:6)%>%filter(TeType!="Unclassified")

ggplot(data_plot, aes(x = depth, y = Percent, color = TeType, shape= taxonomy)) +
  geom_point(alpha=0.6) +
  labs(x = "maximum depth", y = "Percent", color = "TE Type", shape = "Taxonomy") +
  scale_color_manual(values = c("#c91367", "#932598", "#fce2fb", "#fef5ee", "#380a8a")) +
  scale_shape_manual(values = c(17, 16)) +
  theme_classic()+
  geom_abline(intercept=coef(pglsModelDNA)[1], slope= coef(pglsModelDNA)[2], linetype = "dashed", col="#c91367") +
  geom_abline(intercept=coef(pglsModelLINE)[1], slope= coef(pglsModelLINE)[2], linetype = "dashed", col="#932598") +
  geom_abline(intercept=coef(pglsModelLTR)[1], slope= coef(pglsModelLTR)[2], linetype = "dashed", col="#fce2fb") +
  geom_abline(intercept=coef(pglsModelSINE)[1], slope= coef(pglsModelSINE)[2], linetype = "dashed", col="#fef5ee") +
  geom_abline(intercept=coef(pglsModel)[1], slope= coef(pglsModel)[2], linetype = "dashed", col="#380a8a") 

