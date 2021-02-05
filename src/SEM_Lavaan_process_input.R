# This is auxiliary code to load data for the script 
# SEM_Lavaan_NetGroups_2pub.R

# Read Files -------------------------------------------------------------
# --- Paths are automatically set relative to the root of the repo
pathFun=paste(this.dir,"/SEM_analysis/input_data",sep="")
pathInv=pathFun # paste(this.dir,"/main_analysis/input_data",sep="")
pathFunGroup=paste(this.dir,"/functionInk_analysis/output_data",sep="")
pathModel=paste(pathFun,"/LavaanModels/",sep="")
pathOut=paste(this.dir,"/SEM_analysis/output_data",sep="")

#-- Read a single object
if(raw==0){
  setwd(pathFun)
  df.all=readRDS(fileAll)
}else{
  #-- Read in functional data
  setwd(pathFun)
  dd.data=read.csv(fileFun)
  
  #-- Exclude negative controls
  dd.data=dd.data[dd.data$Community!="blank",]
  dd.data=na.omit(dd.data)
  
  #-- Average function measurements across the replicates
  dd=aggregate(dd.data[4:ncol(dd.data)],list(dd.data$Community),function(x){mean(x,na.rm=T)})
  row.names(dd)=unique(dd.data$Community)
  colnames(dd)[1]="Community"
  dd.log=dd
  dd.log[,2:dim(dd)[2]]=log(dd[,2:dim(dd)[2]]+1)
  funTmp=colnames(dd.log)
  Nfun=2*length(funTmp)
  
  #-- Create a vector with the remaining names, select those containing "7" and convert it into a matrix
  idxx=colnames(dd.log) 
  dd.log.sub=subset(dd.log,select=grep("7", idxx))
  dd.log.sub=subset(dd.log,select= -c(pgRPC.7))
  
  #-- Read the  property you prepared to consider for each functional group
  setwd(pathFunGroup)
  funGroups=read.table(fileFunGroup,skip="12",header = TRUE)
  
  #-- Read the invasion file
  setwd(pathInv)
  Invasion.table=read.csv(fileInv)
  inv.log=Invasion.table
  inv.log[,2:dim(inv.log)[2]]=log(inv.log[,2:dim(inv.log)[2]]+1)
  colnames(inv.log)[1]="Community"
  
  # Match communities -------------------------------------------------------
  df.all=merge(dd.log,inv.log) # They both have Community as col name
  #matched=match(dd.log$Community,inv.log$Community)
  df.all=merge(df.all,funGroups,by.x="Community",by.y="Sample")
  #df.all=merge(df.all,Part.dummy,by.x="Community",by.y="Samples")
}



# --- Modify some column names to facilitate model specification

idx=which(colnames(df.all)=="kt2440.cfu.24h") # make life easier changing the colnames
colnames(df.all)[idx]="Putida24h"
idx=which(colnames(df.all)=="kt2440.cfu.96h") # make life easier changing the colnames
colnames(df.all)[idx]="Putida96h"
idx=which(colnames(df.all)=="kt2440.cfu.7d") # make life easier changing the colnames
colnames(df.all)[idx]="Putida168h"
idx=which(colnames(df.all)=="sbw25.cfu.24h") # make life easier changing the colnames
colnames(df.all)[idx]="Fluoresc24h"
idx=which(colnames(df.all)=="sbw25.cfu.96h") # make life easier changing the colnames
colnames(df.all)[idx]="Fluoresc96h"
idx=which(colnames(df.all)=="sbw25.cfu.7d") # make life easier changing the colnames
colnames(df.all)[idx]="Fluoresc168h"

# Rescale data ---------------------------
Nvars=dim(df.all)[2]
all.vars=matrix(0,nrow=1,ncol=Nvars)
sapply(df.all,class)
nofactor=which(sapply(df.all,class)!="factor")
sapply(df.all[,nofactor],var) # check the variances
select.cols=colnames(df.all)[nofactor]
# In ill.scaled covariance matrix, the ratio of the largest to smallest variance 
# is greaterthan say, 100.0 (see e.g. book Rex Kline, pp 81)
# The invasion variables variances are 10-100 larger than the functional groups
# and CPM7 is 10 times larger, so I rescale the vars
rescale10=which(sapply(df.all[select.cols],var)>1)
#rescale100=which(sapply(df.all[select.cols],var)>10)
#colnames(df.all[select.cols][rescale100])# double check 
df.all[select.cols][rescale10]=df.all[select.cols][rescale10]/5 # rescale once 
#df.all[select.cols][rescale100]=df.all[select.cols][rescale100]/5 # a subset rescaled twice
sapply(df.all[,nofactor],var) # check again the variances
sapply(df.all[,nofactor],mean) # check the means