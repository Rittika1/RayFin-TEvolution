# Load the necessary library
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)
library(phytools)
library(geiger)


setwd("/home/rittika/Documents/PolypterusBichir/TEproject/")
# Read data from a file (replace 'your_data_file.csv' with your actual file path)
data <- read.csv('TEdataWithLifeHistory_updatedPCA.csv')
fishtree <- read.tree("/home/rittika/Documents/PolypterusBichir/TEproject/UpdatedTimeTree.phy")
fishtree$tip.label[fishtree$tip.label=="Monocentris_japonicus"] <-"Monocentris_japonica"
fishtree$tip.label[fishtree$tip.label=="Haplochromis_burtoni"] <-"Astatotilapia_burtoni"
fishtree$tip.label[fishtree$tip.label=="Maylandia_zebra"] <-"Metriaclima_zebra"


genomesSize <- data$GenomeSize
names(genomesSize) <- data$Type
TorNOT <- data$teleostOrNot
names(TorNOT) <- data$Type
# name.check(fishtree, genomesSize)
te.ltr <- data$LTR
names(te.ltr) <- data$Type
te.dna <- data$DNA
names(te.dna) <- data$Type
te.line <- data$LINE
names(te.line) <- data$Type
te.sine <- data$SINE
names(te.sine) <- data$Type

# genomesize_totalTE = phylANOVA(fishtree, data$GenomeSize, data$Total_TE)
genomesize_TorNot = phylANOVA(fishtree, TorNOT, genomesSize, posthoc = FALSE)
genomesize_DNA = phylANOVA(fishtree, TorNOT, te.dna, posthoc = FALSE)
genomesize_LTR = phylANOVA(fishtree, TorNOT, te.ltr, posthoc = FALSE)
genomesize_LINE = phylANOVA(fishtree, TorNOT, te.line, posthoc = FALSE)
genomesize_SINE = phylANOVA(fishtree, TorNOT, te.sine, posthoc = FALSE)

##k-test
k.test.genomesize<-phylosig(fishtree, genomesSize, method="K", test=TRUE, nsim=10000)
K.test.DNA<-phylosig(fishtree, te.dna, method="K", test=TRUE, nsim=10000)
K.test.LTR<-phylosig(fishtree, te.ltr, method="K", test=TRUE, nsim=10000)
K.test.LINE<-phylosig(fishtree, te.line, method="K", test=TRUE, nsim=10000)
K.test.SINE<-phylosig(fishtree, te.sine, method="K", test=TRUE, nsim=10000)

#For lambda:
Lambda.test.genomesSize<-phylosig(fishtree, genomesSize, method = "lambda", test=TRUE)
Lambda.test.DNA<-phylosig(fishtree, te.dna, method = "lambda", test=TRUE)
Lambda.test.LTR<-phylosig(fishtree, te.ltr, method = "lambda", test=TRUE)
Lambda.test.LINE<-phylosig(fishtree, te.line, method = "lambda", test=TRUE)
Lambda.test.SINE<-phylosig(fishtree, te.sine, method = "lambda", test=TRUE)

