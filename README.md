# Code for "Phylogeography and paleoclimatic range dynamics explain variable outcomes to contact across a species' range"

This repository provides all code necessary to run down-stream analyses. 


## Major analyses

## Demographic inference

Demographic inference is broken into several sections. These are (1) determining 2N down-sampling, (2) metadata preparation, (3) running Moments, (4) concatenating model runs, and (5) determination of best models and parameter values.

### 2N down-sampling

Down-sampling is done using easySFS (https://github.com/isaacovercast/easySFS). Results should be saved to an Excel file containing the columns (1) POP and (2) BEST, which provides the optimal 2N size.

### metadata preparation

A metadata file is required to run Moments scripts appropriately. To generate, the code file metadata_writer.R is used. It requires files for both the above down-sampled 2N sizes per population as well as the VCF population map (no header).

This metadata file exports a .txt file with columns for (1) iterations - set to 50, (2) Pair - pair name, (3) pop1 - the ID, (4) pop2 - the ID, proj1 - 2N for pop1, proj2 - the ID for pop2.

### running Moments

Numerous model scripts are provided in the /demo_infer/models/ path. For details see (FIG. 2). For models to run appropriately, a VCF (edited directly in *.py scripts) must be provided, as well as the metadata file created above. *.sh files run all models listed in respective *.py files.

These scripts export outputs to a folder (outputs) which must be created prior to running. File names are written as "PAIRNAME_MODEL_expanded_output.txt"

### concatenating models

Models are concatenated in R using the script AIC_concat.R. It binds iterations of all models (e.g., AM, SC...) for all population pairs into a single .csv file (AIC_weights_all_XXperc_linreg_ext*). Note that cyclic models (2C) are concatenated separately from singular models.

### determining the best model and parameters

Files created by concatenating separate runs inform the scripts wAIC.R, which calculates AIC weights and model scores to help judge the best model. Outputs a script named: AIC_weights_all_XXperc.csv

This script combines cyclic and singular model outputs from the concatenate step above, if desired. Outputs of wAIC.R can then be used in the scrpts AICw_Mscore.R and Tsplit_Tmig.R to determine the best runs by AIC weights, AIC, and model scores; as well as plot certain desired parameters (mij/mji; Tmig/Tsplit).


## LEA/sNMF population structure

Accessible in sNMF_structure. Briefly, isolation-by-distance tests should be run (IBD.R) to set expectations about cross-entropy decay and K-clusters. If significant for any lineage, caution should be taken not to inflate K-clusters. 

code in sNMF.R runs the function snmf from LEA on a VCF. Also requires a population map ("noTperf_idpoplin.txt"). this is then saved as a project to be called by sNMF_plot.R. The script can then determine optimal K-clusters in three fashions: (1) quantile slope, (2) DAPC a-score, (3) minimum cross-entropy value.

sNMF_plot.R should be provided a (1) VCF, (2) information from ABBA-BABA tests digested through the script /gene_flow/Dsuite_filter_basicplots.R, (3) group coordinates as a .csv with the columns (a) group - population name, (b) lineage, (c) latitude, (d) longitude. 

A shortened version of this code which can operate given a single Q-matrix can be found in /sNMF_structure/short_code/tess3Q.R. This file only requires a LEA project or a matrix of ancestry coefficients with columns (1) id - individual id, (2) pop, (3) lin - lineage, (4-n) ancestry estimates per cluster. Row totals across these columns should equal 1.


## Other analyses

PCA (/PCA/PCA_full.R) and missing threshold discovery (/VCF_filtering/missing_threhsold.R)

Other scripts are available upon request from the authors.