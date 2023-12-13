library(phytools)
library(geiger)
library(ggplot2)
library(grid)
library(ggrepel)
library("ggphylomorpho")

setwd("/home/rittika/Documents/PolypterusBichir/TEproject/")
fishtree <- read.tree("/home/rittika/Documents/PolypterusBichir/TEproject/UpdatedTimeTree.phy")
fishtree$tip.label[fishtree$tip.label=="Monocentris_japonicus"] <-"Monocentris_japonica"
fishtree$tip.label[fishtree$tip.label=="Haplochromis_burtoni"] <-"Astatotilapia_burtoni"
fishtree$tip.label[fishtree$tip.label=="Maylandia_zebra"] <-"Metriaclima_zebra"
fishtree$tip.label[fishtree$tip.label=="Polypterus_bichir_lapradei"] <-"Polypterus_bichir"
TEtypesdata <- read.csv("/home/rittika/Documents/PolypterusBichir/TEproject/allfishesTEtypes_removed3.csv", header = TRUE, row.names = 1)
TEtypesdata$type <- c("teleost")
TEtypesdata["Polypterus_bichir",]$type <- "non teleost"
TEtypesdata["Polypterus_senegalus",]$type <- "non teleost"
TEtypesdata["Acipenser_ruthenus",]$type <- "non teleost"
TEtypesdata["Lepisosteus_osseus",]$type <- "non teleost"
TEtypesdata["Lepisosteus_oculatus",]$type <- "non teleost"
TEtypesdata["Amia_calva",]$type <- "non teleost"



TEpca <- phyl.pca(fishtree, TEtypesdata, method = "BM", mode = "cov")
biplot(TEpca)
cols <- setNames(c("blue","orange"),levels(TEtypeDNA$type))
phylomorphospace(fishtree, TEtypesdata[,1:2], A=NULL, label=c("radial","horizontal","off"),
                 colors=cols, bty="l",ftype="off",node.by.map=TRUE, node.size=c(0,1.2))


# TEtypeDNA<- TEtypesdata[,c("DNA", "type")]
# rownames(TEtypeDNA) <- rownames(TEtypesdata)
# cols <- setNames(c("blue","orange"),levels(TEtypeDNA$type))
#phylomorphospace(fishtree, TEtypeDNA, A=NULL, label=c("radial","horizontal","off"),
              #   colors=cols, bty="l",ftype="off",node.by.map=TRUE, node.size=c(0,1.2))

# DNApca <- phyl.pca(fishtree, TEtypeDNA, method = "BM", mode = "cov")
# biplot(DNApca)

# ############################
# ##--trying by ggplot
# setdiff(tree$tip.label, matrix$sp)
# ggphylomorpho(fishtree, TEtypesdata, xvar = Dim1, yvar = Dim2,
#                 factorvar =Orden,labelvar = sp)


