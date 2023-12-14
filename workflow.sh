##---These are the steps followed to get the results of the paper: Swimming to Jumping: Unveiling the Dynamic World of Transposable Elements in Ray-finned Fish Evolution.
##-- The file was given in the format .fsa. After an attempt to submit the file to NCBI, the file was rejected because it had many contaminations like adapters and duplicate sequences. 

###---put in the code I used.

### To annotate the genome, I used the following code:  fullbrakerpipeline_bichir.sh
###---RepeatModeler is run as part of the annotation pipeline. This makes two important files: .tbl and .out files. The .tbl file has the percentage of different repeats in the genome. The .out file has the details of the position, length and more details about the repeats. 
##-- getTEdetails.py takes in .tbl file and print TE details in the tabular format

##-- TEcontentVSgenomesize.py take in the tbl file and the genome size and gives the percentage of repeats in the genome.

##-- getTEdetails_everyfamily.py takes in the .out file and prints the details of each family in the tabular format.

##---Pfam analysis was done using the following code: Pfam-analysis.sh. The results were plotted using R script.
## boxplot.R    
## LinearModelsForTEs.R       
## phylopca_bichir_AD.R  -- makes the phylomorphospace
## TE_types_allfishesplot.R   -- calculates pgls and makes the plot
## t-testandlambatest_bichir.R -- calculates Pagels lambda and Blombergs K
## LinearModels_for_TEs_OU.R  -- calculates OU models
## model_fitting_TEclasses.R  -- phylostep model for stepwise pgls model
## TreePrAbHeatmap.R -- SIMMAP code for ancestral state reconstruction