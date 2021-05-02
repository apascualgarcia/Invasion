#!/usr/bin/perl -w
# *******************************************
# * SamplePropsOFtaxaClusters.pl *
# *******************************************
#
# This script aims to characterize abundances and number of species, when in the data there is
# a partition in taxa (we call them clusters) , obtained for instance
# with NodeLinkage.pl (functionInk) after computing their correlations.
# For each sample, it identifies the taxa that are found in specific taxa clusters, 
# and compute properties of the taxa for the whole cluster, such as the maximum abundance
# of the species in the cluster, maximum relative abundance, median of the z.score, etc. These properties will
# be associated to that specific cluster, for instance the total abundance of taxa in cluster i in sample a. 
# The idea is to reduce the complexity of the species in the different samples talking about clusters of taxa. 
# This is the script I use to perform next structural equation modelling.
# NOTE: This script is an extension of Samples2TaxaClusters.pl which, in turn, is a modification of  TaxaClusters2function.pl.
# INPUTS: The file with the clusters and a list of the cluster indices you want to process, that should be filled in the code below. It follows the format of NodeLinkage.pl
#         A matrix of taxons (columns) Vs samples (rows) containing raw abundances (not resampled)
# OUTPUTS: One file with different properties of the taxa clusters in the samples, one line per sample:
#                 $AvSp,$StdvSp,$AvSpFrac,$StdvSpFrac,$AvAb,$StdvAb,$AvFracAb,$StdvFracAb
#          One file with the fraction of abundances for each taxa cluster
#########################################
# Theoretical biology (ETH)
# April 16th, 2020. Alberto Pascual-García 
# alberto.pascual.garcia@gmail.com 
#########################################

# -- Additional modules
use POSIX; # Used for floor and ceil functions
use Statistics::Basic qw(:all); # mean, median, etc, opt1
use PDL; # compute mean, etc, opt2
use List::Util qw( min max );
   
# -- Determine the clusters you want to work with


@clusters=("1","2","3","5","8","13","14","16","17","18","20","22","23","24","27","28","35","39","45"); # SparCC AvLink Step 91 with Paenibacillus species merge

@labels=("Yellow","Red","Green","Orange","LightGreen","OliveBrown","Pink","Black"); # Identify clusters with colors in conseunsus set 


# -- Input paths. Here paths are relative to the output directory, from which I run the script
$pathClust='../functional_groups_detection/'; # SparCC all (functionInk)  SimilarityMetaComm/'; # sparCC consensus between metacommunities
$pathAbund='../../../SparCC/input_data/';

# -- Input files
$fileClust='Clusters-NL_Average_StopStep-91_Network_BasisCor_SparCC_GT0.2_pval0.01_FunGroups_Time0_Paenibacillus.txt';# File with the clusters
$fileAbund='composition_otu_matched.csv';# Note that we normalize below so by number of reads

# -- Output files

$Set='Time0';
$TaxaClus='NL_Average_StopStep-91'; # You may want to consider the labels of the input files here
$genericOutput='SamplePropsOFtaxaClus_'.$Set.'_'.$TaxaClus.'_';
$pathOut=''; #$pathClust;

# -- Starting computation
print "  \n";
print "*********************************************************  \n";
print "* Quantifying species properties for different clusters *  \n";
print "*********************************************************  \n";
print "  \n";

# -- Open output files and build up headers
&filesOutSetUp();

# -- Read the taxa found in each cluster
print "~~~ Read the taxa found in each cluster  \n";

$fileIn=$pathClust.$fileClust;
open(INTMP,$fileIn)  or die "Can't find file $fileIn\n";
@INTMP = <INTMP>;
$i=0; # This index will control the clusters
foreach $Clus (@clusters){ # new version, make independent of order April 2020
    foreach $line(@INTMP){
	chomp($line);
	#$Clus=$clusters[$i]; # Clusters are in increasing order in the file --> not any more, April 2020
	$string="CLUS_".$Clus."_NODES"; # Build the identifier for the first cluster
	@fields=split(" ",$line);
	if($fields[0] eq $string){ # look for the linens where the nodes are defined (the taxa)
	    print  join(' ','.. Processing cluster/',$Clus,' in line: ',$string),"\n";
	    foreach $Taxa (@fields){ # Process the taxa one by one
		if($Taxa eq $string){next;} # The first element is still the identifier
		push(@{$Clus2taxa{$Clus}},$Taxa);
		#$Taxa2clus{$Taxa}=$Clus; # CHANGE WITH PREVIOUS VERSION: Taxa is only in one cluster, 
		#$ClustersList{$Clus}=1
	    }
	    $i+=1;
	    #if($i == $#clusters+1){
		#last;
	    #}
	}
    }
}


# -- Look for the samples in which each taxa has been found (Abundances matrix).
print "~~~ Looking in which samples each cluster taxa was found  \n";
$fileIn=$pathAbund.$fileAbund;
open(INTMP,$fileIn)  or die "Can't find file $fileIn\n";
@INTMP = <INTMP>;
$i=0; # This index will control that we are in the first line
foreach $line (@INTMP){
    chomp($line);
    @fields=split(",",$line);
    if($i == 0){ # In the first line we find the name of the N taxa
	$i+=1;
	$j=0;
	foreach $Taxa(@fields){
	    if($Taxa eq "Community"){next;} # skip the first name
	    $j+=1; # We start storing in one because in the next lines there are N+1 columns 
	    $Taxa2index[$j]=$Taxa; # Vector of strings
	    $CumAbund{$Taxa}=0;
	    $CumAbund2{$Taxa}=0;
	    $TaxaList{$Taxa}=1;
	    #print  join(' ','.. Processing taxa:',$Taxa),"\n";
	}
    }else{
	$j=0;
	$SpTmp=0;
	$AbTmp=0;
	foreach $obs (@fields){ # Now we identify the sample and then check which taxa were observed
	    if($j==0){
		$Sample=$obs; # The first column is the name of the sample
		$SamplesList{$Sample}=1;
	    }else{
		$obs=log($obs+1); # Transform the abundances
		$Taxa=$Taxa2index[$j]; # Identify the taxa associated to that index
		$CumAbund{$Taxa}+=$obs;
		$CumAbund2{$Taxa}+=$obs**2;
		$Abundances{$Sample}{$Taxa}=$obs; # REMEMBER Log-transform abundances!
		if($obs != 0){ # it means that in Sample the $j taxa was observed ($obs>0)
		    #push(@{$Taxa2Sample{$Taxa}},$Sample); # Associate the sample to the taxa
		    #push(@{$Sample2Taxa{$Sample}},$Taxa); # Associate the taxa to the sample	   
		    $SpTmp+=1;
		    $AbTmp+=$obs;
		    #print  join(' ','.. Taxa ',$Taxa,' observed in sample ',$Sample),"\n"; # DEBUG
		}		  
	    }
	    $j+=1;
	} # -- end foreach observation
	$TotalSpecies{$Sample}=$SpTmp;
	$TotalAbundance{$Sample}=$AbTmp;
    }
}
@Samples=keys%SamplesList;
@TaxaKeys=keys%TaxaList;
print "~~~~ I've found $#Samples+1 Samples and $#TaxaKeys+1 Taxa  \n";

# -- Compute averages and stdv of the log-transformed taxa abundances for Zscores
print "~~~ Compute average and standard deviation abundances  \n";
foreach $Taxa (@TaxaKeys){
    $MeanAbund{$Taxa}=$CumAbund{$Taxa}/($#Samples+1);
    $StdvAbund{$Taxa}=sqrt($CumAbund2{$Taxa}/($#Samples+1)-$MeanAbund{$Taxa}**2);
}

# -- Assign species and abundances to every sample/cluster and samples to cluster/partition
print "~~~ Compute clusters properties  \n";

# ---- We need first to build a header for the output file
$u=-1;
foreach $metric (@metrics){ # The metrics are defined in te function filesOutSetUp
    $u+=1;
    $ASA=$Handles[$u];
    print $ASA 'Sample',"\t";
    foreach $Cluster (@clusters){
	print $ASA 'C'.$Cluster,"\t";
    }
    print $ASA "Total \n";
}

foreach $Sample (@Samples){
    $u=-1;
    foreach $metric (@metrics){ # Print the name of the sample in each file first
	$u+=1;
	$ASA=$Handles[$u];	
	print $ASA $Sample,"\t";
    }
    $AbundAllClus=0;
    foreach $Cluster (@clusters){
	$SpeciesTmp=0;
	$AbundTmp=0;
	@ZscoreVec = ();
	@AbundVec = ();
	@RelAbVec = ();
	foreach $Taxa (@{$Clus2taxa{$Cluster}}){ 
	    #print join(" ",'> DEBUG',$Sample,$Taxa,$Abundances{$Sample}{$Taxa},$MeanAbund{$Taxa},$StdvAbund{$Taxa}),"\n"; # DEBUG
	    $ZscoreTmp=($Abundances{$Sample}{$Taxa}-$MeanAbund{$Taxa})/$StdvAbund{$Taxa};
	    $AbundTmp=$Abundances{$Sample}{$Taxa};
	    $RelAbTmp=$AbundTmp/$TotalAbundance{$Sample};
	    push(@ZscoreVec,$ZscoreTmp); # store values of taxa within each cluster
	    push(@AbundVec,$AbundTmp); 
	    push(@RelAbVec,$RelAbTmp); 
	    $SpeciesCumTmp+=1;	       
	    $AbundCumTmp+=$AbundTmp;
	}
	$AbundAllClus+=$AbundTmp;
	$ClusterMetric{$metrics[0]}=$SpeciesCumTmp;
	$ClusterMetric{$metrics[1]}=$AbundCumTmp;
	my $piddle = pdl @AbundVec;
	($mean,$prms,$median,$min,$max,$adev,$rms) = statsover $piddle; # extract several metrics
	$ClusterMetric{$metrics[2]}=$max; # (@AbundVec);
	$ClusterMetric{$metrics[3]}=$mean; # (@AbundVec);
	$ClusterMetric{$metrics[4]}=$median; # (@AbundVec);
	my $piddle = pdl @RelAbVec;
	($mean,$prms,$median,$min,$max,$adev,$rms) = statsover $piddle; # extract several metrics
	$ClusterMetric{$metrics[5]}=$max; # (@RelAbVec);
	$ClusterMetric{$metrics[6]}=$mean; # (@RelAbVec);
	$ClusterMetric{$metrics[7]}=$median; # (@RelAbVec);
	my $piddle = pdl @ZscoreVec;
	($mean,$prms,$median,$min,$max,$adev,$rms) = statsover $piddle; # extract several metrics
	$ClusterMetric{$metrics[8]}=$max; # (@ZscoreVec);
	$ClusterMetric{$metrics[9]}=$mean; # (@ZscoreVec);
	$ClusterMetric{$metrics[10]}=$median; # (@ZscoreVec);
	# print join(" /",'> DEBUG',@ZscoreVec,$ClusterMetric{$metrics[10]}),"\n"; # DEBUG

	if($SpeciesTmp==0){$SpeciesTmp=1;} # this is just to avoid zero division, no further consequences
	$u=-1;
	 foreach $metric (@metrics){
	     $u+=1;
	     $ASA=$Handles[$u];
	     print $ASA $ClusterMetric{$metric},"\t";
	 }

    }
    $Remainder=$TotalAbundance{$Sample}-$AbundAllClus; # These are the abundances not attributable to clusters
    $TotalAbundMean=mean(values %TotalAbundance);
    #print join(" /",'> DEBUG',values %TotalAbundance),"\n"; # DEBUG
    #print join(" /",'> DEBUG',$TotalAbundMean),"\n"; # DEBUG
    #exit; # DEBUG

    $TotalAbundStdv=stddev(values %TotalAbundance);
    $ZscoreTotal=($TotalAbundance{$Sample}-$TotalAbundMean)/$TotalAbundStdv;
    #print join(" /",'> DEBUG',$Sample,$TotalSpecies{$Sample},$TotalAbundMean,$TotalAbundStdv,$TotalAbundance{$Sample}, $ZcoreTotal),"\n"; # DEBUG

    $u=-1;
    foreach $metric (@metrics){
	$u+=1;
	$ASA=$Handles[$u];
	if($u==0){
	    $TotalOut=$TotalSpecies{$Sample};
	}elsif($u<4){
	    $TotalOut=$TotalAbundance{$Sample};
	}elsif($u<7){
	    $TotalOut=1;
	}else{
	    $TotalOut=$ZscoreTotal;
	}
	print $ASA $TotalOut,"\n";    
    }
}


print '~~~ Done!',"\n";
print '  ',"\n";
print '****************************',"\n";
print '** Program finished',"\n";
print '** Check your results  ',"\n";
print '** Bye! ',"\n";
print '****************************',"\n";
print '  ',"\n";
print '  ',"\n";




###################################################
##                                       FUNCTIONS
###################################################


######################
#        filesOutSetUp
######################
# Define handlers for output files, open them and
# creates a header with information of parameters and date. It
# returns the array of handlers.

sub filesOutSetUp{
    $Nfiles=2; # Files needed 
    &theTime();
    @metrics=("SpCount","AbundCum","AbundMax","AbundMean","AbundMedian","RelAbMax","RelAbMean","RelAbMedian","ZscoreMax","ZscoreMean","ZscoreMedian");
    $u=-1;
    foreach $metric (@metrics){
	$u+=1;
	$ASA=*OUT.$u; #
	$particularOutput=$metric.".dat";
	push(@Handles,$ASA);	    
	$file=$pathOut.$genericOutput.$particularOutput;
	print '>> Opening file: ',$file,"\n";
	open($ASA, ">$file") || die "Couldn't open file $file"; 
	print $ASA '# >> Output from SampleProps2TaxaClusters.pl <<',"\n";
	print $ASA '# >> Input file for the list of clusters: ',$fileClust,"\n";
	print $ASA '# >> Input file for samples to taxa matrix: ',$fileAbund,"\n";
	print $ASA '# >> Identifiers of the clusters: ',"\n";
	print $ASA join(' ','# >> : ',@clusters),"\n";
	print $ASA '# >> Labels of the clusters: ',"\n";
	print $ASA join(' ','# >> : ',@labels),"\n";
	print $ASA '# >> This output file: ',$file,"\n";
	print $ASA '# >> Running at date: ',$Time,"\n";
	print $ASA '# >> Theoretical biology, ETH-Zürich, A.P-G.', "\n";
	print $ASA '#WARNING: Abundances values are LOG-transformed ', "\n";
	print $ASA '# ', "\n";
	#print $ASA '#';	
    }
    return @Handles,@metrics;
}

######################
#     theTime
######################
# Return the time where the script is being executed

sub theTime{
    ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset) = localtime();
    $year = 1900 + $yearOffset;
    $Time = " (mm,dd,yy) $month, $dayOfMonth, $year, and time: $hour:$minute:$second,";
    return $Time;
}
