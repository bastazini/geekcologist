# Darwin to the Rescue: Using Phylogenetic information to overcome the Raunkiaeran shortfall
#teste
###Load packages
require(geiger)
require(PVR)
require(phytools)
require(missForest)  
library(ape)
library(phytools)
library(picante)

####Import trait.df data from github
##data comes from Vizentin‐Bugoni, Jeferson, et al. "Including rewiring in the estimation of the robustness of mutualistic networks." Methods in Ecology and Evolution 11.1 (2020): 106-116.
#make sure you select "raw"option in github, before copiyng the webpage direction
trait.df= read.csv("https://raw.githubusercontent.com/vanderleidebastiani/rewiring/master/DataSetExamples/SantaVirginia/SantaVirginia_dataset_h_morph.csv",
                 header = TRUE,sep=",",row.names=1)


###Here we will create missing values,
#you can also use the function of missforest package for that
#missForest::prodNA()

###Define the proportion of missing observations 
miss.val=3

while(sum(is.na(trait.df) == TRUE) < (miss.val)){
  trait.df[sample(nrow(trait.df),1), sample(ncol(trait.df),1)] = NA
}
#Look at the new dataset with missing values
trait.df


###Import phylogentic information (nexus file)  with 100 phylogenies from  GitHub directory
#Tree distribution from: "The global diversity of birds in space and time; W. Jetz, G. H. Thomas, J. B. Joy, K. Hartmann, A. O. Mooers, doi:10.1038/nature11631]
#Subsampled and pruned from birdtree.org on 2021-11-03

tree=read.nexus("https://raw.githubusercontent.com/bastazini/geekcologist/main/birds_phylo.nex")

##Create a consensus tree. 
#Note: p=0.5 specifies that the tree must be "majority rules consensus (MRC)"

hb.tree=consensus.edges(tree, consensus.tree=consensus(tree,p=0.5))  

consensus_ultra=chronos(hb.tree, lambda=0)  
# Lambda is the rate-smoothing parameter.  Generally, we want lambda to be small, so there is very little smoothing  See here (citation at bottom): http://www.justinbagley.org/1226/update-new-functions-for-generating-starting-trees-for-beast-or-starbeast-in-r

##Plot phylogenetic tree
plotTree(hb.tree,type="fan")
  


###Imputting trait.df data using the phylogenetic information
# this is based on the  example from github (https://github.com/vanderleidebastiani/missForestImputation)

## Imputation
# Decompose the phylogenetic distance matrix into a set of orthogonal vectors (PVRs)
  
phylo.vectors = PVR::PVRdecomp(hb.tree)

# Extract the PVRs
pvrs = phylo.vectors@Eigen$vectors

# Combine trait.dfs and PVRs in data frame
trait.dfs.pvrs = cbind(trait.df, pvrs)

# Imputation using missForest (note that this function have other arguments, see details in Stekhoven and Buehlmann 2012)
RF.imp = missForest::missForest(trait.dfs.pvrs, maxiter = 15, ntree = 100, variablewise = FALSE)


# Here it is! Your imputed dataset!
trait.dfs.imp = RF.imp$ximp[, seq_len(ncol(trait.df)), drop = FALSE]
trait.dfs.imp

## Acess imputation error using  Normalized Root-Mean-Square Deviation (NRMSD)
# NRMSE ranges from 0 to ≈ 1. 
# NRMSE values  ≈ 1 occur when the estimations are poor or when the noise involved is too large
RF.imp$OOBerror


