library(phytools)
library(geiger)
library(ggplot2)
library(dplyr)

#only load libraries you are actually using, I think the above are all you need at this stage (might be wrong, I still have R open from the workshop on friday)
#library(grid) 
#library(ggrepel)

#I don't know what you are doing here.... 
#library("ggphylomorpho")

#Before we even start, here is a best practice. Please define all your paths as objects here. It is my pet peeve to have to read an entire script to figure out what files I need and also looks sloppy. I will all caps next time, I have nightmares of debugging 1000s of lines of code from people with dozens of files that appear with no warning.... 

#paths
wd<-"/home/rittika/Documents/PolypterusBichir/TEproject/"
tree_path<-"/home/rittika/Documents/PolypterusBichir/TEproject/NoNodeLabelsUpdatedTimeTree.phy"
data_path<-"/home/rittika/Documents/PolypterusBichir/TEproject/allfishesTEtypes_removed3.csv"
TEDetails_data_path<-"/home/rittika/Documents/PolypterusBichir/TEproject/TEDetails_only4family.csv"



#with paths defined, we can just focus on the code below and once it is finalized it does not need to be edited :-)

#set working directory
setwd(wd)

#read in data/tree
fishtree <- read.tree(tree_path)
TEtypesdata <- read.csv(data_path, header = TRUE, row.names = 1)

#modify tip labels from timetree.org
fishtree$tip.label[fishtree$tip.label=="Monocentris_japonicus"] <-"Monocentris_japonica"
fishtree$tip.label[fishtree$tip.label=="Haplochromis_burtoni"] <-"Astatotilapia_burtoni"
fishtree$tip.label[fishtree$tip.label=="Maylandia_zebra"] <-"Metriaclima_zebra"
fishtree$tip.label[fishtree$tip.label=="Polypterus_bichir_lapradei"] <-"Polypterus_bichir"

#the data has total TE as a column, obviously this is the sum of the other columns so we don't want to include it as a variable. We also do not want the type at this stage since we won't use it until we plot much later
PCAdata<-TEtypesdata[,1:5]

#we want to make sure our tree and data match before going further. This is important, do not skip these steps. 

## CHECK THAT THE NAMES ALL MATCH BETWEEN TREE AND DATA, this drops bathygadus from the data
td<-treedata(fishtree, PCAdata) #this function is in the geiger package

updatedPCAdata <-as.data.frame(td$data)

#We can now do a phylogenetic PCA of the %s
TEpca <- phyl.pca(fishtree, updatedPCAdata, method = "BM", mode = "cov")
#get the % variation per axis
summary(TEpca)

#now get the actual PC scores for each species.IMPORTANT. THIS IS WHAT YOU NEED TO PLOT. Your original code just plotted the raw data!
TEScores<-TEpca$S

#Read in this phylomorphospace to ggplot convert function 
#phylomorpho to dataframe function, modified from phytools blog
phytools2ggplot<-function(phylomorphospace){
  output <-data.frame(
    xstart= phylomorphospace$xx[phylomorphospace$edge[,1]],
    ystart=phylomorphospace$yy[phylomorphospace$edge[,1]],
    xstop=phylomorphospace$xx[phylomorphospace$edge[,2]],
    ystop=phylomorphospace$yy[phylomorphospace$edge[,2]],
    nodestart=phylomorphospace$edge[,1],
    nodestop=phylomorphospace$edge[,2])
  return(output)
}

#Draw the phylomorphospace, save as a ggplot-able object in the process
TESpace<-phytools2ggplot(phylomorphospace(fishtree, TEScores[,1:2],label="off",xlab="PCA1",ylab="PCA2") )

#Now we can add convex hulls. We will start by defining two colors here                 
colors <- c("orange", "blue")

#define groups, I added erpetoichthys, this was missing again!
plotscores<-as.data.frame(TEScores)
plotscores$type <- c("teleost")
plotscores["Polypterus_bichir_lapradei",]$type <- "nonteleosts"
plotscores["Polypterus_senegalus",]$type <- "nonteleosts"
plotscores["Erpetoichthys_calabaricus",]$type <- "nonteleosts"
plotscores["Acipenser_ruthenus",]$type <- "nonteleosts"
plotscores["Lepisosteus_osseus",]$type <- "nonteleosts"
plotscores["Lepisosteus_oculatus",]$type <- "nonteleosts"
plotscores["Amia_calva",]$type <- "nonteleosts"

p <- ggplot(plotscores, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors)
  #view it
  p
  #now that we have a normal scatterplot, we can draw the convex hulls
 for (i in 1:length(unique(plotscores$type))) {
  hull_points <- plotscores[plotscores$type == sort(unique(plotscores$type))[i], 1:2]
  hull <- chull(hull_points)
  hull_points <- hull_points[c(hull, hull[1]), ]
  p <- p + geom_polygon(data = hull_points, aes(x = PC1, y = PC2), fill=colors[i],color = colors[i], alpha = 0.2)
}
#view it
p
#now we can add the tree using geom_segment since it is just lines in a space
p+ geom_segment(data = TESpace, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.4, color = "gray") + theme_classic()

#############################################
#Part 2 repeat with the  zoomed in Data
#############################################

#read in data
#Part 2.A PCA of the entire matrix

All_data<-read.csv(TEDetails_data_path)[-1]
All_data_Names<-read.csv(TEDetails_data_path)[1]
rownames(All_data)<-All_data_Names[,1]
All_data<-All_data %>% select(-c(DNA, LTR, LINE,SINE))

tdad<-treedata(fishtree, All_data) #this function is in the geiger package

updatedPCAdata.2 <- as.data.frame(tdad$data) 
updatedPCAdata.2[is.na(updatedPCAdata.2)] <- 0
#We can now do a phylogenetic PCA of the %s
TEpca.2 <- phyl.pca(tdad$phy, updatedPCAdata.2, method = "BM", mode = "cov")
#get the % variation per axis
summary(TEpca.2 )

#now get the actual PC scores for each species.IMPORTANT. THIS IS WHAT YOU NEED TO PLOT. Your original code just plotted the raw data!
TEScores.2<-TEpca.2$S

#Draw the phylomorphospace, save as a ggplot-able object in the process
TESpace.2<-phytools2ggplot(phylomorphospace(tdad$phy, TEScores.2[,1:2],label="off",xlab="PCA1",ylab="PCA2") )

#Now we can add convex hulls. We will start by defining two colors here                 
colors <- c("orange", "blue")

#define groups, I added erpetoichthys, this was missing again!
plotscores<-as.data.frame(TEScores.2)
plotscores$type <- c("teleost")
plotscores["Polypterus_bichir_lapradei",]$type <- "nonteleosts"
plotscores["Polypterus_senegalus",]$type <- "nonteleosts"
plotscores["Erpetoichthys_calabaricus",]$type <- "nonteleosts"
plotscores["Acipenser_ruthenus",]$type <- "nonteleosts"
plotscores["Lepisosteus_osseus",]$type <- "nonteleosts"
plotscores["Lepisosteus_oculatus",]$type <- "nonteleosts"
plotscores["Amia_calva",]$type <- "nonteleosts"

pdf("phylomorpho-allTE.pdf")
p <- ggplot(plotscores, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors)
  #view it
  p
  #now that we have a normal scatterplot, we can draw the convex hulls
 for (i in 1:length(unique(plotscores$type))) {
  hull_points <- plotscores[plotscores$type == sort(unique(plotscores$type))[i], 1:2]
  hull <- chull(hull_points)
  hull_points <- hull_points[c(hull, hull[1]), ]
  p <- p + geom_polygon(data = hull_points, aes(x = PC1, y = PC2), fill=colors[i],color = colors[i], alpha = 0.2)
}
#view it
p
#now we can add the tree using geom_segment since it is just lines in a space
p+ geom_segment(data = TESpace.2, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.4, color = "gray") + theme_classic()

dev.off()
#################################
###Repeat with just DNA
################################
updatedPCAdata.DNA<-updatedPCAdata.2 %>% select(contains("DNA"))

#We can now do a phylogenetic PCA of the %s
TEpca.DNA <- phyl.pca(tdad$phy, updatedPCAdata.DNA, method = "BM", mode = "cov")
#get the % variation per axis
summary(TEpca.DNA)

#now get the actual PC scores for each species.
TEScores.DNA<-TEpca.DNA$S

#Draw the phylomorphospace, save as a ggplot-able object in the process
TESpace.DNA<-phytools2ggplot(phylomorphospace(tdad$phy, TEScores.DNA[,1:2],label="off",xlab="PCA1-DNA",ylab="PCA2-DNA") )

#Now we can add convex hulls. We will start by defining two colors here                 
colors <- c("orange", "blue")

#define groups, I added erpetoichthys, this was missing again!
plotscores<-as.data.frame(TEScores.DNA)
plotscores$type <- c("teleost")
plotscores["Polypterus_bichir_lapradei",]$type <- "nonteleosts"
plotscores["Polypterus_senegalus",]$type <- "nonteleosts"
plotscores["Erpetoichthys_calabaricus",]$type <- "nonteleosts"
plotscores["Acipenser_ruthenus",]$type <- "nonteleosts"
plotscores["Lepisosteus_osseus",]$type <- "nonteleosts"
plotscores["Lepisosteus_oculatus",]$type <- "nonteleosts"
plotscores["Amia_calva",]$type <- "nonteleosts"

p <- ggplot(plotscores, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors)
  #view it
  p
  #now that we have a normal scatterplot, we can draw the convex hulls
 for (i in 1:length(unique(plotscores$type))) {
  hull_points <- plotscores[plotscores$type == sort(unique(plotscores$type))[i], 1:2]
  hull <- chull(hull_points)
  hull_points <- hull_points[c(hull, hull[1]), ]
  p <- p + geom_polygon(data = hull_points, aes(x = PC1, y = PC2), fill=colors[i],color = colors[i], alpha = 0.2) 
}
#view it
p
#now we can add the tree using geom_segment since it is just lines in a space
p+ geom_segment(data = TESpace.DNA, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.4, color = "gray") + theme_classic() 
                 
#################################
###Repeat with just LINE
################################
updatedPCAdata.LINE<-updatedPCAdata.2 %>% select(contains("LINE"))

#We can now do a phylogenetic PCA of the %s
TEpca.LINE <- phyl.pca(tdad$phy, updatedPCAdata.LINE, method = "BM", mode = "cov")
#get the % variation per axis
summary(TEpca.LINE)

#now get the actual PC scores for each species.
TEScores.LINE<-TEpca.LINE$S

#Draw the phylomorphospace, save as a ggplot-able object in the process
TESpace.LINE<-phytools2ggplot(phylomorphospace(tdad$phy, TEScores.LINE[,1:2],label="off",xlab="PCA1",ylab="PCA2") )

#Now we can add convex hulls. We will start by defining two colors here                 
colors <- c("orange", "blue")

#define groups, I added erpetoichthys, this was missing again!
plotscores<-as.data.frame(TEScores.LINE)
plotscores$type <- c("teleost")
plotscores["Polypterus_bichir_lapradei",]$type <- "nonteleosts"
plotscores["Polypterus_senegalus",]$type <- "nonteleosts"
plotscores["Erpetoichthys_calabaricus",]$type <- "nonteleosts"
plotscores["Acipenser_ruthenus",]$type <- "nonteleosts"
plotscores["Lepisosteus_osseus",]$type <- "nonteleosts"
plotscores["Lepisosteus_oculatus",]$type <- "nonteleosts"
plotscores["Amia_calva",]$type <- "nonteleosts"

p <- ggplot(plotscores, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors)
  #view it
  p
  #now that we have a normal scatterplot, we can draw the convex hulls
 for (i in 1:length(unique(plotscores$type))) {
  hull_points <- plotscores[plotscores$type == sort(unique(plotscores$type))[i], 1:2]
  hull <- chull(hull_points)
  hull_points <- hull_points[c(hull, hull[1]), ]
  p <- p + geom_polygon(data = hull_points, aes(x = PC1, y = PC2), fill=colors[i],color = colors[i], alpha = 0.2)
}
#view it
p
#now we can add the tree using geom_segment since it is just lines in a space
p+ geom_segment(data = TESpace.LINE, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.4, color = "gray") + theme_classic()                 

#################################
###Repeat with LTR
################################
updatedPCAdata.LTR<-updatedPCAdata.2 %>% select(contains("LTR"))

#We can now do a phylogenetic PCA of the %s
TEpca.LTR <- phyl.pca(tdad$phy, updatedPCAdata.LTR, method = "BM", mode = "cov")
#get the % variation per axis
summary(TEpca.LTR)

#now get the actual PC scores for each species.
TEScores.LTR<-TEpca.LTR$S

#Draw the phylomorphospace, save as a ggplot-able object in the process
TESpace.LTR<-phytools2ggplot(phylomorphospace(tdad$phy, TEScores.LTR[,1:2],label="off",xlab="PCA1",ylab="PCA2") )

#Now we can add convex hulls. We will start by defining two colors here                 
colors <- c("orange", "blue")

#define groups, I added erpetoichthys, this was missing again!
plotscores<-as.data.frame(TEScores.LTR)
plotscores$type <- c("teleost")
plotscores["Polypterus_bichir_lapradei",]$type <- "nonteleosts"
plotscores["Polypterus_senegalus",]$type <- "nonteleosts"
plotscores["Erpetoichthys_calabaricus",]$type <- "nonteleosts"
plotscores["Acipenser_ruthenus",]$type <- "nonteleosts"
plotscores["Lepisosteus_osseus",]$type <- "nonteleosts"
plotscores["Lepisosteus_oculatus",]$type <- "nonteleosts"
plotscores["Amia_calva",]$type <- "nonteleosts"

p <- ggplot(plotscores, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors)
#view it
p
#now that we have a normal scatterplot, we can draw the convex hulls
for (i in 1:length(unique(plotscores$type))) {
  hull_points <- plotscores[plotscores$type == sort(unique(plotscores$type))[i], 1:2]
  hull <- chull(hull_points)
  hull_points <- hull_points[c(hull, hull[1]), ]
  p <- p + geom_polygon(data = hull_points, aes(x = PC1, y = PC2), fill=colors[i],color = colors[i], alpha = 0.2)
}
#view it
p
#now we can add the tree using geom_segment since it is just lines in a space
p+ geom_segment(data = TESpace.LTR, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                linewidth = 0.4, color = "gray") + theme_classic()    


#################################
###Repeat with just SINE
################################
updatedPCAdata.SINE<-updatedPCAdata.2 %>% select(contains("SINE"))

#We can now do a phylogenetic PCA of the %s
TEpca.SINE <- phyl.pca(tdad$phy, updatedPCAdata.SINE, method = "BM", mode = "cov")
#get the % variation per axis
summary(TEpca.SINE)

#now get the actual PC scores for each species.
TEScores.SINE<-TEpca.SINE$S

#Draw the phylomorphospace, save as a ggplot-able object in the process
TESpace.SINE<-phytools2ggplot(phylomorphospace(tdad$phy, TEScores.SINE[,1:2],label="off",xlab="PCA1",ylab="PCA2") )

#Now we can add convex hulls. We will start by defining two colors here                 
colors <- c("orange", "blue")

#define groups, I added erpetoichthys, this was missing again!
plotscores<-as.data.frame(TEScores.SINE)
plotscores$type <- c("teleost")
plotscores["Polypterus_bichir_lapradei",]$type <- "nonteleosts"
plotscores["Polypterus_senegalus",]$type <- "nonteleosts"
plotscores["Erpetoichthys_calabaricus",]$type <- "nonteleosts"
plotscores["Acipenser_ruthenus",]$type <- "nonteleosts"
plotscores["Lepisosteus_osseus",]$type <- "nonteleosts"
plotscores["Lepisosteus_oculatus",]$type <- "nonteleosts"
#plotscores["Amia_calva",]$type <- "nonteleosts"

p <- ggplot(plotscores, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors)
#view it
p
#now that we have a normal scatterplot, we can draw the convex hulls
for (i in 1:length(unique(plotscores$type))) {
  hull_points <- plotscores[plotscores$type == sort(unique(plotscores$type))[i], 1:2]
  hull <- chull(hull_points)
  hull_points <- hull_points[c(hull, hull[1]), ]
  p <- p + geom_polygon(data = hull_points, aes(x = PC1, y = PC2), fill=colors[i],color = colors[i], alpha = 0.2)
}
#view it
p
#now we can add the tree using geom_segment since it is just lines in a space
p+ geom_segment(data = TESpace.SINE, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                linewidth = 0.4, color = "gray") + theme_classic()    

# 
# ##################################
# ###Repeat with just longnose gar
# ################################
#All_data<-read.csv(TEDetails_data_path)[-1]
All_data_Names<-read.csv(TEDetails_data_path)[1]
rownames(All_data)<-All_data_Names[,1]
All_data<-All_data %>% select(-c(DNA, LTR, LINE,SINE))

tipstodrop <- c("Amia_calva", "Polypterus_bichir", "Polypterus_senegalus" , "Lepisosteus_oculatus", "Acipenser ruthenus" , "Erpetoichthys_calabaricus")

fishtree <- drop.tip(fishtree, tipstodrop)
tdad<-treedata(fishtree, All_data) #this function is in the geiger package

updatedPCAdata.5 <- as.data.frame(tdad$data) 
updatedPCAdata.5[is.na(updatedPCAdata.5)] <- 0
#We can now do a phylogenetic PCA of the %s
TEpca.5 <- phyl.pca(tdad$phy, updatedPCAdata.5, method = "BM", mode = "cov")
#get the % variation per axis
summary(TEpca.5 )

#now get the actual PC scores for each species.IMPORTANT. THIS IS WHAT YOU NEED TO PLOT. Your original code just plotted the raw data!
TEScores.5<-TEpca.5$S

#Draw the phylomorphospace, save as a ggplot-able object in the process
TESpace.5<-phytools2ggplot(phylomorphospace(tdad$phy, TEScores.5[,1:2],label="off",xlab="PCA1",ylab="PCA2") )

#Now we can add convex hulls. We will start by defining two colors here                 
colors <- c("orange", "blue")

#define groups, I added erpetoichthys, this was missing again!
plotscores<-as.data.frame(TEScores.5)
plotscores$type <- c("teleost")
# plotscores["Polypterus_bichir_lapradei",]$type <- "nonteleosts"
# plotscores["Polypterus_senegalus",]$type <- "nonteleosts"
# plotscores["Erpetoichthys_calabaricus",]$type <- "nonteleosts"
# plotscores["Acipenser_ruthenus",]$type <- "nonteleosts"
plotscores["Lepisosteus_osseus",]$type <- "nonteleosts"
# plotscores["Lepisosteus_oculatus",]$type <- "nonteleosts"
# plotscores["Amia_calva",]$type <- "nonteleosts"

p <- ggplot(plotscores, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors)
#view it
p
#now that we have a normal scatterplot, we can draw the convex hulls
for (i in 1:length(unique(plotscores$type))) {
  hull_points <- plotscores[plotscores$type == sort(unique(plotscores$type))[i], 1:2]
  hull <- chull(hull_points)
  hull_points <- hull_points[c(hull, hull[1]), ]
  p <- p + geom_polygon(data = hull_points, aes(x = PC1, y = PC2), fill=colors[i],color = colors[i], alpha = 0.5)
}
#view it
p
#now we can add the tree using geom_segment since it is just lines in a space
p+ geom_segment(data = TESpace.5, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                linewidth = 0.4, color = "gray") + theme_classic()

