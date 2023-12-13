#triggerfish raincloud plots

##load libraries
library(tidyverse)
library(ggridges)
library(viridis)

#set file paths
triggerDataPath<-"~/Documents/PolypterusBichir/TEproject/TEdataWithLifeHistory_updatedPCA.csv"

#read in data
data<-read.csv(triggerDataPath)

#write a plotting function to keep the code lines down
###############################################
###Read in this function, do not modify 
###############################################
# habitatcloud<-function(data, adultTrait="Adult_Dorsal_aspect_ratio", juvTrait="Juvi_Dorsal_aspect_ratio", plotName="Dorsal Aspect Ratio",...){
# 	#wrangle data for ggplot
# 	targetAdultData<-data %>% select(c(adultTrait, "X3.adult.habitat")) %>% mutate(Lifestage="adult") %>% rename(Trait= adultTrait, Habitat=X3.adult.habitat)
# 	targetjuvData<-data %>% select(c(juvTrait,"X3.juvi.habitat")) %>% mutate(Lifestage="juvenile")  %>% rename(Trait= juvTrait, Habitat=X3.juvi.habitat)
# 	#combine
# 	ggplotData<-rbind(targetAdultData, targetjuvData) %>% na.omit() %>% unite(combinedFactor, Lifestage:Habitat)
# 	#plot the raincloud plots
# 	ggplot(ggplotData, aes(x= Trait, y= combinedFactor, fill = ..x..))+ 
# 	geom_density_ridges_gradient(scale = 0.85, rel_min_height = 0.02, jittered_points = TRUE, position="raincloud",quantile_lines=TRUE, quantile_fun=function(Trait,...)mean(Trait)) +
#   scale_fill_viridis(name = plotName, option = "C", direction = 1) + 
#   theme_classic()
# }

habitatcloud<-function(data, adultTrait="teleostOrNot", plotName="TEplot name random",...){
  #wrangle data for ggplot
  targetAdultData<-data %>% select(c(adultTrait, "teleostOrNot")) %>% mutate(Lifestage="teleost") %>% rename(Trait= adultTrait, Habitat=X3.adult.habitat)
  targetjuvData<-data %>% select(c(juvTrait,"teleostOrNot")) %>% mutate(Lifestage="NonTeleost")  %>% rename(Trait= juvTrait, Habitat=X3.juvi.habitat)
  #combine
  ggplotData<-rbind(targetAdultData, targetjuvData) %>% na.omit() %>% unite(combinedFactor, Lifestage:Habitat)
  #plot the raincloud plots
  ggplot(ggplotData, aes(x= Trait, y= combinedFactor, fill = ..x..))+ 
    geom_density_ridges_gradient(scale = 0.85, rel_min_height = 0.02, jittered_points = TRUE, position="raincloud",quantile_lines=TRUE, quantile_fun=function(Trait,...)mean(Trait)) +
    scale_fill_viridis(name = plotName, option = "C", direction = 1) + 
    theme_classic()
}


###############################################
#ok, modify here for each set of traits following the first example
###############################################

#Dorsal Incidence angle
habitatcloud(data, adultTrait="Adult_dorsal_incidence_angle", juvTrait="Juvi_dorsal_incidence_angle", plotName="Dorsal Incidence Angle")

habitatcloud(data, adultTrait="Adult_Body_RW1", juvTrait="Juvi_Body_RW1", plotName="Body RW 1")

habitatcloud(data, adultTrait="Adult.A.J.Dorsal1", juvTrait="Juv.A.J.Dorsal1", plotName="Dorsal RW 1")

