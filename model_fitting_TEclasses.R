library(phylolm)


##------------EXAMPLE CODE---------
# set.seed(123456)
# tre = rcoal(60)
# taxa = sort(tre$tip.label)
# b0=0; b1=1;
# x1 = rTrait(phy=tre,model="BM",
#             parameters=list(ancestral.state=0,sigma2=10))
# x2 = rTrait(phy=tre,model="BM",
#             parameters=list(ancestral.state=0,sigma2=10))
# x3 = rTrait(phy=tre,model="BM",
#             parameters=list(ancestral.state=0,sigma2=10))
# y <- b0 + b1*x1 + 
#   rTrait(n=1,phy=tre,model="BM",parameters=list(
#     ancestral.state=0,sigma2=1))
# dat = data.frame(trait=y[taxa],pred1=x1[taxa],pred2=x2[taxa],pred3=x3[taxa])
# fit = phylostep(trait~pred1+pred2+pred3,data=dat,phy=tre,model="BM",direction="both")
# summary(fit)

##--define filepaths here
wd<-"/home/rittika/Documents/PolypterusBichir/TEproject/"
TEdata_path<-"/home/rittika/Documents/PolypterusBichir/TEproject/TEdataWithLifeHistory_updatedPCA.csv"
fishtree_Path<-"~/Documents/PolypterusBichir/TEproject/NoNodeLabelsUpdatedTimeTree.phy"

##--load in the TE file
TEtypesdata <- read.csv(TEdata_path, header = TRUE, row.names=1)
##--subset only the four types of TEs
TEclasses <- TEtypesdata[, 1:4]

##-- read in tree
fishtree <- read.tree(fishtree_Path)

##--calculate fit using phylostep
phylostep(DNA~LTR + teleostOrNot +Habitat, data = TEtypesdata, phy = fishtree, model = "BM", direction = "both")
phylostep(DNA~LINE+ teleostOrNot +Habitat, data = TEtypesdata, phy = fishtree, model = "BM", direction = "both")
phylostep(DNA~SINE + teleostOrNot +Habitat, data = TEtypesdata, phy = fishtree, model = "BM", direction = "both")
phylostep(LTR~LINE + teleostOrNot +Habitat, data = TEtypesdata, phy = fishtree, model = "BM", direction = "both")
phylostep(LTR~SINE + teleostOrNot +Habitat, data = TEtypesdata, phy = fishtree, model = "BM", direction = "both")
phylostep(LINE~SINE + teleostOrNot +Habitat, data = TEtypesdata, phy = fishtree, model = "BM", direction = "both")
###---- teleost or not vs LTR+LINE+SINE+DNA
# y = DNA
# x = LTR + taxonomy

