README
=======


## Description

* author = Alberto Pascual-García
* web = apascualgarcia.github.io
* date = June 24th, 2020 (ETH-Zürich)
* description = This script takes the invasion data [1], a processed matrix with information regarding the abundances of functional groups found in [2] for the experiments described in [3, 4], and the productivities of the communities, and it explores with structural equation models the relation between composition (encoded in the functional groups), productivity and invasion.
* edit = The user must define the SEM models in lavaan format in external files located in the folder lavaanModels.Then, the variable selectModel should point to that file. The input data should be located in the respective directories. The paths to the file relative to the repository are provided, paths should be edited if data is stored in another location. Note that the function to retrieve the files works only with *Rstudio*, and it should be manually included otherwise.
* usage =If only a single model is run, the model should be launched from the script: `launch_SEM_single.R`. If multiple models are run, the model should be launched from the script: `launch_SEM_multiple.R`. See details in the respective scripts.
* results = The script will run each model and store it in a .mod file, and a summary of results of the fit and the MI indexes is stored in another file (.fit). It will then create a representation with the script `SEM_PathPlot.R` for the plain model (excluding variances), the structural model and the full model. If several models are analysed with the "multiple" version, an analysis of all the models can be made.

## REFERENCES

[1] Jones, M. L., Rivett, D. W., Pascual-García, A., & Bell, T. (2019). Productive bacterial communities exclude invaders. bioRxiv.

[2] Pascual‐García, A., & Bell, T. (2020). functionInk: An efficient method to detect functional groups in multidimensional networks reveals the hidden structure of ecological communities. Methods in Ecology and Evolution, 11(7), 804-817.

[3] Pascual-García, A., & Bell, T. (2020). Community-level signatures of ecological succession in natural bacterial communities. Nature communications, 11(1), 1-11.

[4] Rivett, D. W. & Bell, T. Abundance determines the functional role of bacterial phylotypes in complex communities. Nat. Microbiol. 3, 767 (2018).


