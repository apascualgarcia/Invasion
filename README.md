README
=======


## Description

* author = Alberto Pascual-García
* web = apascualgarcia.github.io
* date = June 24th, 2020 (ETH-Zürich)
* description = This repository contains a pipeline to estimate the influence of microbial composition and community-level function for the experiments described in [3,4] in invasion resistance, presented in [1].


## Pipeline

All the scripts used in the pipeline are in the directory `src`. Input and output of data for each step are within a directory with the name of the step.

1. Estimation of bacterial co-occurrences and segregations

  * **script**: `network_inference_sparCC.R`
  * **input**: Working directory `SparCC/input_files`.
   * A pre-processed otu matrix neglecting rare OTUs (see [1] for details): `composition_otu_matched.csv`
  * **output**: Working directory `SparCC/output_files`.
   * Basis correlations estimated by SparCC (wide format): `BasisCor-wide_SparCC_FunGroups_Time0.txt`.
   * Basis correlations and pseudo pvalues in long format: `BasisCor-Pval-long_SparCC_FunGroups_Time0.txt`
   * Final network contaning only the correlations that passed the filter with the parameters (correlation greater than `GT$value` and pvalue lower than the values indicated in the label): `Network_BasisCor_SparCC_GT$val_pval$val_FunGroups_Time0.txt`

2. Determination of functional groups

  * **script**: Community detection performed with functionink [2], available in [this repository](https://github.com/apascualgarcia/functionInk)
  * **input**: The network obtained in the previous step.
  * **output**: In directory `functionInk/output_data/functional_groups_detection/`
   * All the output provided by functionink when the clustering is run with no stopping criteria "no stop"
   * A figure showing the partition densities, in which it is shown that the maximum of the total partition density is achieved at step 91. 
   * The output of functionink at the step 91, in which we retrieve the clusters. 
   * A Cytoscape project with the network, after filtering out functional groups with less than three members but including _Paenibacillus_ species `Network_BasisCor_SparCC_GT0.2_pval0.01_FunGroups_Time0.cys`. A png file with the final network is also included.
   * A file with the functional groups selected `selected_clusters_$param.txt` and README file explaining the criteria to select the functional groups.

3. Computation of values

  * **script**: `SamplePropsOFtaxaClusters.pl`
  * **input**: The OTU table, located in `SparCC/input_files` and the file `Clusters-NL_Average_StopStep-91_Network_BasisCor_SparCC_GT0.2_pval0.01_FunGroups_Time0_Paenibacillus.txt` in `functionInk_analysis/output_data/functional_groups_detection/`,
  * **output**: In directory `functionInk/output_data/functional_groups_values/`, different files with different metrics for the functional groups.

4. Structural equation modelling

  * **scripts**:
   * edit = The user must define the SEM models in lavaan format in external files located in the folder `lavaanModels`.Then, the variable `selectModel` should point to that file. The input data should be located in the respective directories. The paths to the file relative to the repository are provided, paths should be edited if data is stored in another location. Note that the function to retrieve the files works only with *Rstudio*, and it should be manually included otherwise.
   * usage =If only a single model is run, the model should be launched from the script: `launch_SEM_single.R`. If multiple models are run, the model should be launched from the script: `launch_SEM_multiple.R`. Other functions are coded in different files, see details in the respective scripts.
  * **input**: 
   * The output from previous step: `functionInk/output_data/functional_groups_values/SamplePropsOFtaxaClus_Time0_NL_Average_StopStep-91_ZscoreMean.dat`
   * The invasion experiments: `SEM_analysis/input_data/invasionexp_jan21.csv`
   * The functional measurements: `SEM_analysis/input_data//20151016_Functions_remainder.csv`
   * Once the analysis is run for the first time, the script creates a RDS file with all the previous inputs processed, and can be used this file reducing the computational burden: `SEM_analysis/input_data/semdf.RDS`.
  * **output**: `SEM_analysis/output_data/`
    * The script will run each model and store it in a .mod file, and a summary of results of the fit and the MI indexes is stored in another file (.fit). It will then create a representation with the script `SEM_PathPlot.R` for the plain model (excluding variances), the structural model and the full model. If several models are analysed with the "multiple" version, an analysis of all the models can be made.

## REFERENCES

[1] Jones, M. L., Rivett, D. W., Pascual-García, A., & Bell, T. (2019). Productive bacterial communities exclude invaders. bioRxiv.

[2] Pascual‐García, A., & Bell, T. (2020). functionInk: An efficient method to detect functional groups in multidimensional networks reveals the hidden structure of ecological communities. Methods in Ecology and Evolution, 11(7), 804-817.

[3] Pascual-García, A., & Bell, T. (2020). Community-level signatures of ecological succession in natural bacterial communities. Nature communications, 11(1), 1-11.

[4] Rivett, D. W. & Bell, T. Abundance determines the functional role of bacterial phylotypes in complex communities. Nat. Microbiol. 3, 767 (2018).


