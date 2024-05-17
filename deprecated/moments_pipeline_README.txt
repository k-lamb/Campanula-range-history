Moments pipeline explained

##########################################################################################################

For simple models:

Simulations run on Rivanna
Processing run on own system

Rivanna section:

main scripts directory
/home/...

to run each model, see sub-folders, then run: All_moments_{model}.py script with the shell run_dadi_pairs_{MODEL}.sh
outputs go into: ./All/{model}/outputs/

To concatenate:
first run: AIC_concat.sh which grabs all results and compiles from each model subfolder using AIC_moments.py and AIC_metadata.txt
It grabs lowest AIC run for each model, translates raw parameters, then concatenates in a single file: 
./All/{model}/outputs/{model}_aic_max_best_translated.txt

then it will also spit out a file with the best run of model (for a given model) of each pair of populations in a single file:
./All/{model}/{model}_best_sorted.csv

Home system section:

download csv's produced for each model and place in: /Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/
add column with model name in excel to each
sort by name of pair then add A, W, or B as the 'group' column manually
then run: /Users/kericlamb/Desktop/Documents/Research/Paper Code/Debban-Lamb/Moments/moments/model_concat.py
R analysis done in: /Users/kericlamb/Desktop/Documents/Research/Paper Code/Debban-Lamb/Moments/moments/graphs/graphs.R

######################################################################################################

For extended models (from Momigliano et al.):

to run each model, use All_moments_{model}_expanded.py and edit respective *.sh file for use (see MIG1 for an example)

run extended_AIC_concat.R with file directories in setwd sprintf clause... 

these models contain the simple models (MIG, SC, NO_MIG), so no need to run simple model structure shown above unless you want one_pop models (NULL)

To calculate weights, use AIC_weights_graphs.R