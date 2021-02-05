#################################################
#      SEM models with Lavaan                   #
#################################################
# author = Alberto Pascual-García
# web = apascualgarcia.github.io
# date = June 24th, 2020 (ETH-Zürich)
# description = This script is a wrapper of SEM_Lavaan_NetGroups_2pub.R and
#     it permits to run it several times with different models. However, if you
#     are interested in comparing several models, since the anova
#     function below to analyse different SEM has a rigid sintax, the script cannot be used without recoding the call of
#     this function. 
########## 
#rm(list=ls())

#-- Load packages needed
library(gplots) 
library(MASS)
require(lavaan)
library(rstudioapi)    

#### START EDITING
# Models and files --------------------------------------------------------
# --- Provide the files for your models in a list
# Basic model
#models=c("CompProdInvasion_NoMediation.lav",
#         "CompProdInvasion_PartialMediation.lav","CompProdInvasion_CompleteMediation.lav")
# Basic model with diversity
models=c("CompProdInvasion_NoMediation_Diversity.lav",
         "CompProdInvasion_PartialMediation_Diversity.lav",
         "CompProdInvasion_CompleteMediation_Diversity.lav")
# Optimized model
#models=c("CompProdInvasion_NoMediation_mod6.2.10.lav",
#         "CompProdInvasion_PartialMediation_mod6.2.10.lav",
#         "CompProdInvasion_CompleteMediation_mod6.2.10.lav")
# --- Determine if input data should be loaded (required first time the code is run)
raw=1 # = 1 if process the data from scracth, 0 read a pre-processed RDS file

#--- Input files for raw = 1
# .... The functions and invasion experiment
fileFun="20151016_Functions_remainder.csv"
fileInv="invasiondata_albertomatched_050617.csv"

# .... The properties of the functional groups you will use for your model 
fileFunGroup="SamplePropsOFtaxaClus_Time0_NL_Average_StopStep-90_ZscoreMean.dat"

# .... Alternatively, load a single R data object (raw = 0)
fileAll="semdf.RDS"

# .... The output files saved with the model are created according to "selectModel" variable)

# ... Comparison of models
comparison=1 # if you want to compare models (=1), =0 otherwise. Note that this part of the code should be recoded (see below)
###### STOP EDITING

model.fit.list=list()
k=0
for(selectModel in models){
  k=k+1
  cat("** Running model ",selectModel,"\n")
  if(k == 1){
    this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] 
    pathCode=paste(this.dir,"/src",sep="") # this the path where the code is
    setwd(pathCode)
    source("SEM_Lavaan_process_input.R")
  }
  setwd(pathCode)
  source("SEM_Lavaan_NetGroups_2pub.R")
  model.fit.list[[k]]=fit
}

if(comparison == 1){
  # Here you should choose the models you want to compare
  # ..... For nested models
  # 1="Independent",2="NoMediation",3="Partial",4="Complete"
  # df(Partial) < df(Complete) = df(Independent) < df(NoMed)
  anova.test=lavTestLRT(model.fit.list[[2]],model.fit.list[[3]]) # Partial and complete
  
  
  # ..... Rank by parsimony (AIC), we expect a deltaAIC > 2 per degree of freedom lost
  library(AICcmodavg) 
  aic.rank=aictab(model.fit.list,
                  modnames=c("NoMediation","Partial","Complete"))
  
}

setwd(pathOut)
sink("MultipleTests_results.txt")
cat("******** \n")
cat("Anova tests \n")
cat("1=NoMediation,2=Partial,3=Complete \n")
print(anova.test)
cat(" \n")
cat("******** \n")
cat("Rank by parsimony (AIC), we expect a deltaAIC > 2 per degree of freedom lost \n")
print(aic.rank)
sink()
