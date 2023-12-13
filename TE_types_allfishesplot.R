library(ggplot2)
library(ggdendro)
library(ape)
library(geiger)
library(nlme)
library(phytools)

TEvsGenomedata <- read.csv("/home/rittika/Documents/PolypterusBichir/TEproject/TEvsGenome_updatedwithAvatree.csv", header = TRUE, row.names = 1)
TEtypesdata <- read.csv("/home/rittika/Documents/PolypterusBichir/TEproject/allfishesTEtypes_removed3.csv", header = TRUE, row.names = 1)
TEtypesdata$GenomeLength <- TEvsGenomedata$Total_Length[match(rownames(TEtypesdata), rownames(TEvsGenomedata))]
TEtypesdata$type <- c("teleost")
TEtypesdata["Polypterus_bichir",]$type <- "non teleost"
TEtypesdata["Polypterus_senegalus",]$type <- "non teleost"
TEtypesdata["Acipenser_ruthenus",]$type <- "non teleost"
TEtypesdata["Lepisosteus_osseus",]$type <- "non teleost"
TEtypesdata["Lepisosteus_oculatus",]$type <- "non teleost"
TEtypesdata["Amia_calva",]$type <- "non teleost"

TEtypesdata$GenomeLength <- TEtypesdata$GenomeLength/1000000
##--write this to a file for safekeeping later
write.csv(TEtypesdata, "/home/rittika/Documents/PolypterusBichir/TEproject/allfishdetails.csv", row.names = TRUE)

##-- scatter plot with the TE content percent and the genome size data
plot(log(TEtypesdata$GenomeLength), TEtypesdata$Total_TE, col = factor(TEtypesdata$type), pch=19)
barplot(TEtypesdata$GenomeLength, col = factor(TEtypesdata$type))

##--read in the tree and modify the tip labels
fishtree <- read.tree("/home/rittika/Documents/PolypterusBichir/TEproject/UpdatedTimeTree.phy")
fishtree$tip.label[fishtree$tip.label=="Monocentris_japonicus"] <-"Monocentris_japonica"
fishtree$tip.label[fishtree$tip.label=="Haplochromis_burtoni"] <-"Astatotilapia_burtoni"
fishtree$tip.label[fishtree$tip.label=="Maylandia_zebra"] <-"Metriaclima_zebra"

##-- modify data to run pgls
TEpercent <- TEtypesdata$Total_TE
Genomesize <- TEtypesdata$GenomeLength

##--give them names
names(TEpercent) <- names(Genomesize) <- rownames(TEtypesdata)

##--pgls with total TE
pglsModel <- gls(Total_TE ~ GenomeLength, correlation = corBrownian(1,fishtree),
                 data = TEtypesdata, method = "ML")
summary(pglsModel)
coef(pglsModel)
plot(Total_TE~GenomeLength, data=TEtypesdata, col = factor(TEtypesdata$type), pch=19)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])

par(mfrow=c(2,2))
##--pgls with DNA
pglsModelDNA <- gls(DNA ~ GenomeLength, correlation = corBrownian(1,fishtree),
                 data = TEtypesdata, method = "ML")
summary(pglsModelDNA)
coef(pglsModelDNA)
plot(DNA~GenomeLength, data=TEtypesdata, col = factor(TEtypesdata$type), pch=19, main = "DNA ~ Genomesize")
abline(a = coef(pglsModelDNA)[1], b = coef(pglsModelDNA)[2])

##--pgls with LINE
pglsModelLINE <- gls(LINE ~ GenomeLength, correlation = corBrownian(1,fishtree),
                    data = TEtypesdata, method = "ML")
summary(pglsModelLINE)
coef(pglsModelLINE)
plot(LINE~GenomeLength, data=TEtypesdata, col = factor(TEtypesdata$type), pch=19, main = "LINE ~ Genomesize")
abline(a = coef(pglsModelLINE)[1], b = coef(pglsModelLINE)[2])

##pgls with SINE
pglsModelSINE <- gls(SINE ~ GenomeLength, correlation = corBrownian(1,fishtree),
                    data = TEtypesdata, method = "ML")
summary(pglsModelSINE)
coef(pglsModelSINE)
plot(SINE~GenomeLength, data=TEtypesdata, col = factor(TEtypesdata$type), pch=19, main = "SINE ~ Genomesize")
abline(a = coef(pglsModelSINE)[1], b = coef(pglsModelSINE)[2])

##--pgls with LTR
pglsModelLTR <- gls(LTR ~ GenomeLength, correlation = corBrownian(1,fishtree),
                    data = TEtypesdata, method = "ML")
summary(pglsModelLTR)
coef(pglsModelLTR)
plot(LTR~GenomeLength, data=TEtypesdata, col = factor(TEtypesdata$type), pch=19, main = "LTR ~ Genomesize")
abline(a = coef(pglsModelLTR)[1], b = coef(pglsModelLTR)[2])

# 
# for ( i in colnames(TE))
