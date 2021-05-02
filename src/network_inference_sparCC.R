##################################################
# network_inference.R
##################################################
# In this script I infer correlations between OTUs or samples abundances
# with sparCC, I recompute the same metric with randomized networks,
# and finally print some networks with significant links and corr > abs(0.3).
# -------
# Note: This script is computationally intensive both from cpu and
# memory usage. As a reference, for a matrix 600x600, 4 bootstrap
# realizations take ~50 mins in 2 cpus. This is was, however, in a test
# with a laptop with 8GB and 4 cpus only. A possibly relevant increase
# in the time may come from the computer swapping after running short in RAM.
# It  is better to be conservative in the use of resources even if it takes
# longer, because if a cpu fails at returning results the script will halt with
# no results.
# ......
# Lausanne, April 24th, 2021
# Theoretical Biology, ETH
# apascualgarcia.github.io
###################################################

library(SpiecEasi)
library(phyloseq)
library(reshape2)

###### START EDITING
# ... directories for input  and output data
dirInRel="SparCC/input_data" # directory for the otu table relative to the root dir of the repo
dirOutRel="SparCC/output_data" # directory for output networks relative to the root dir of the repo
fileOTU="composition_otu_matched.csv" # name of the otu table
# ... computational params
nboot=200 #999 # number of bootstraps used to address the confidence interval
ncpus=10 # number of cpus used in the parallelization
# ... outputs
thr=0.2 # threshold to determine that correlation is meaningful
pval0=0.01 # pvalue to determine that a correlation is significant
label="FunGroups_Time0" # an optional label to identify your data
###### STOP EDITING

# --- Set directories
#this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] # root  dir of the repo
                # note that if you run this script in a cluster the above command only works in rstudio
this.dir="/home/apascual/Nextcloud/Research/Projects/FunctionalGroups/Repositories/Invasion"
dirIn=paste(this.dir,dirInRel,sep="/")
dirOut=paste(this.dir,dirOutRel,sep="/")

# --- Load otu table
setwd(dirIn)
otu.in=read.csv(fileOTU,row.names = 1)#,skip=1)
head(otu.in)[1:5,1:5] # double check
#colnames(otu.in)=t(column.names[2:length(column.names)]) # observe transposition here
#otu.in=t(otu.in) # to compute correlations between samples we transpose here

# Start the computation --------
# --- Compute the basis correlations
sparcc=sparcc(otu.in) # Compute the correlations and covariances
sparcc.cor=sparcc$Cor # get the correlations
colnames(sparcc.cor)=colnames(otu.in)
rownames(sparcc.cor)=colnames(otu.in)
#sparcc.cor=as.matrix(read.table(file=fileOut)) # used after a correct computation but bad formatting

# --- Bootstrap the correlations
# ..... test running time
# system.time({sparccBoot=sparccboot(otu.in,R=4,ncpus = 2)}) # 8GB RAM
# user   system  elapsed 
#3202.796   11.174 2377.700 
sparccBoot=sparccboot(otu.in,R=nboot,ncpus = ncpus)
sparccPval=pval.sparccboot(sparccBoot) # retrieve empirical pvalues from bootstrapped distro
sparccPval.df=as.data.frame(sparccPval)

# Write outputs --------
setwd(dirOut)
fileOut=paste("BasisCor-wide_SparCC_",label,".txt",sep="")
write.table(sparcc.cor,file=fileOut,quote = FALSE)

# --- Transform into network format
diag(sparcc.cor) = -11 # Identify diagonal and lower triangle to remove them below
sparcc.cor[lower.tri(sparcc.cor)]= -11
sparcc.graph=melt(sparcc.cor) # transform into long matrix, only upper triangle needed
sparcc.graph=sparcc.graph[which(sparcc.graph$value != -11),]
# sparcc.graph=data.frame() # transform manually into long matrix
# for(i in 1:(dim(sparcc.cor)[1]-1)){
#   for(j in (i+1):dim(sparcc.cor)[1]){
#     df.tmp=c(colnames(sparcc.cor)[i],colnames(sparcc.cor)[j],
#                         sparcc.cor[i,j])
#     df.tmp=as.matrix(df.tmp)
#     df.tmp=as.data.frame(t(df.tmp))
#     sparcc.graph=rbind(sparcc.graph,df.tmp)
#   }
# }
colnames(sparcc.graph)=c("OTU1","OTU2","corr")

# .... filter those with a low correlation
id.sig=which((abs(sparccPval.df$cors) > thr) & (sparccPval.df$pvals < pval0))
sparccPval.sig.df=sparccPval.df[id.sig,]
sparcc.graph.sig=sparcc.graph[id.sig,] # get just those significant
#plot(sparcc.graph.sig$corr,sparccPval.sig.df$cors) # double check both df are the same

# .... filter those non-significant
sparcc.graph.sig$type=matrix(data=NA,nrow=dim(sparcc.graph.sig)[1],ncol=1)
sparcc.graph.sig$type[which(sparcc.graph.sig$corr>0)]=1
sparcc.graph.sig$type[which(sparcc.graph.sig$corr<0)]=-1
fileOut=paste("Network_BasisCor_SparCC_GT",thr,"_pval",pval0,"_",label,".txt",sep="")  
write.table(sparcc.graph.sig,file=fileOut,quote=FALSE,row.names = FALSE,sep="\t")

# .... print the pvalues
sparcc.pval=data.frame(sparcc.graph[,1:2],sparccPval.df)
colnames(sparcc.pval)=c("OTU1","OTU2","BasisCor","Pval")
fileOut=paste("BasisCor-Pval-long_SparCC_",label,".txt",sep="")  
write.table(sparcc.pval,file=fileOut,quote=FALSE,row.names = FALSE,sep="\t")
