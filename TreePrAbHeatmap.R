#load libraries
library("tidyverse")
library("ggtree")
library("ggtreeExtra")
library("caper")
library("parallel")
library("hrbrthemes")

#file paths
wd<-"~/Documents/TE_Evolution/ASR_Parts"
tree_path<-"NoNodeLabelsUpdatedTimeTree.phy"
data_path<-"presence-absencedata.csv"

#set wd
setwd(wd)

#read in data
teleostTree <- ggtree::read.tree(tree_path)
teRawData <- read.csv(data_path, row.names=1, header = TRUE)

#add names to track data
rownames(teRawData)<-teRawData[,1]
justvalues<-teRawData[,-1]

#td test
test<-justvalues[,1]
names(test)<-rownames(teRawData)

td<-treedata(teleostTree, test)

p<-ggtree(td$phy) + 
    geom_tiplab(size=2, align=TRUE, linesize=.5) + 
    theme_tree2()
heatMapPAM(p, justvalues, offset=100, width=0.6, 
        colnames=TRUE, colnames_angle=45, colnames_offset_y = 2, font.size = 2, colnames_position = "top", col_colours="blue") +
    scale_x_ggtree() + 
    scale_y_continuous(expand=c(0, 0.3))

#########################################
#### SIMMAP ASR
########################################

#AIC functions for model fitting from phytools blog
#first one calculates AIC score
aic.phytools<-function(logL,k) 2*k-2*logL

#this one calculates aic weight
aic.w.phytools<-function(aic){
    d.aic<-aic-min(aic)
    exp(-1/2*d.aic)/sum(exp(-1/2*d.aic))
}

Rootgenome<-matrix(nrow=25, ncol=2)
HolostTeleogenome<-matrix(nrow=25, ncol=2)
JustTeleogenome<-matrix(nrow=25, ncol=2)


for (i in 2:26){
print(i)
x<-teRawData[,i] 
names(x)<-rownames(teRawData)
rowtarget<-i-1
#hard coding an exception, circle back to make this general
if(i ==14){
remove<-which(is.na(x)==TRUE)
x<-x[-remove]
td<-treedata(teleostTree, x) }

#if all present we don't need to do an ASR/can't
if(sum(x)==107){
Rootgenome[rowtarget,]<-c(0,1)
HolostTeleogenome[rowtarget,]<-c(0,1)
JustTeleogenome[rowtarget,]<-c(0,1)
} else {

#Find the best model for transitions in state
logL<-sapply(c("ER","SYM","ARD"),
             function(model,tree,x) make.simmap(tree,x,model)$logL,
             tree=td$phy,x=x)
#AIC summaries
AIC<-mapply(aic.phytools,logL,c(1,3,6))
AIC
AIC.W<-aic.w.phytools(AIC)
AIC.W

#From this we can normalize to the total number of simulations we will do to incorporate model uncertainty
nsim<-5000
Nsim<-round(nsim*AIC.W)
d<-if(sum(Nsim)>nsim) -1 else 1
nsim<-Nsim+d*sample(c(rep(1,abs(nsim-sum(Nsim))),
                      rep(0,length(Nsim)-abs(nsim-sum(Nsim)))))
nsim

nsim<-nsim[nsim!=0]
trees<-NULL
for(j in 1:length(nsim)){
  obj<-make.simmap(td$phy,x,model=names(nsim)[j],nsim=nsim[j])
  if(nsim[j]==1){ 
    obj<-list(obj)
    class(obj)<-"multiPhylo"
  }
  trees<-c(trees,obj)
}
class(trees)<-"multiPhylo"
pd<-describe.simmap(trees)
Rootgenome[rowtarget,]<-pd$ace[1,]
HolostTeleogenome[rowtarget,]<-pd$ace[4,]
JustTeleogenome[rowtarget,]<-pd$ace[5,]
print(pd$ace[1,])
}
}   

#plot the three ancestral mobilomes!
TEs<-colnames(teRawData)[-1]
data<-tibble(TEs= TEs, Actinopterygii= Rootgenome[,2], Neopterygii= HolostTeleogenome[,2], Teleostei= JustTeleogenome[,2])
pivotData<-data %>% pivot_longer(cols=2:4)

ggplot(pivotData, aes(TEs, name, fill= value)) + 
  geom_tile(color="white", lwd=1, linetype=1) +
  scale_fill_gradient(low="white", high="blue") +
  theme_ipsum() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
