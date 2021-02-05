#################################################
#     launch single SEM model with Lavaan       #
#################################################
# author = Alberto Pascual-García
# web = apascualgarcia.github.io
# date = June 24th, 2020 (ETH-Zürich)
# description = This code launches a single SEM analysis
#  specified by the user. 
rm(list=ls())

#-- Load packages needed
library(gplots) 
library(MASS)
require(lavaan)
library(rstudioapi)    

#### START EDITING
# Models and files --------------------------------------------------------
# --- Provide the file for your model
selectModel="CompProdInvasion_Independent.lav"

# --- Determine if input data should be loaded (required first time the code is run)
data.required=1 # =1 reading input data is needed = 0, data is already loaded.
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
###### STOP EDITING

# read the input data if required
if(data.required == 1){
  this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] 
  pathCode=paste(this.dir,"/src",sep="") # this the path where the code is
  setwd(pathCode)
  source("SEM_Lavaan_process_input.R")
}

# run the model
source("SEM_Lavaan_NetGroups_2pub.R")

