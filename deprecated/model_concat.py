import pandas as pd
import numpy as np
from natsort import natsort_keygen

mig = pd.read_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/MIG1_best_sorted.csv')
sc = pd.read_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/SC1_best_sorted.csv')
no_mig = pd.read_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/NO_MIG_best_sorted.csv')
am = pd.read_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/AM_best_sorted.csv')
one_pop = pd.read_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/ONE_POP_best_sorted.csv')

#only do this if working from new one_pop file. this removes a duplication issue
one_pop = one_pop.drop(columns=['index'])
one_pop = one_pop.drop_duplicates(subset= 'pair',keep='first', ignore_index=True)
one_pop.to_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/ONE_POP_best_sorted.csv')

#appends all df's into one df
df = sc.append(mig, sort=False)
df = df.append(no_mig, sort=False)
df = df.append(one_pop, sort=False)
df = df.append(am, sort=False)
df = df.drop(columns=['index']) #drop an old index column 

df.to_csv("/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/All_models_results.csv")

############################################################################################################################################
#
# AIC model selection. Important! There are two ways, best model by AIC and best model by AIC +/- 2pts
#
############################################################################################################################################

# #keep only best model
best = df.sort_values('aic').drop_duplicates(subset= 'pair', keep='first')
best.to_csv("/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/BEST_models_results.csv")

#keep best model and all within 2 AIC points use script below:
pairs = pd.read_csv("/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/All_pairs.txt", header = None) #reads in pairs list
pairs = pairs[0].values.tolist() #converts df to list

#set up df that will be read into 
PMmod=open('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/best_AIC.csv','w')
PMmod.write(
    str("index")+','+
    str("pair")+','+
    str("nu1")+','+
    str("nu2")+','+
    str("Ts")+','+
    str("Tsc")+','+
    str("Tm")+','+
    str("m12")+','+
    str("m21")+','+
    str("theta")+','+
    str("ll_model")+','+
    str("AIC")+','+
    str("nu1_num")+','+
    str("nu2_num")+','+
    str("mji")+','+
    str("mij")+','+
    str("mig_pop1")+','+
    str("mig_pop2")+','+
    str("dv_time")+','+
    str("sc_time")+','+
    str("model")+'\n'
    )
PMmod.close()

#subsets each pair and grabs best model by AIC and all within 2 AIC points of it, then writes to df created above.
for i in range(len(pairs)):
    subset = df[df["pair"] == pairs[i]] #subset to each pair
    subset = subset[subset.aic>(subset.aic.min()+2)] #grab models within 2 AIC 
    subset.to_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/best_AIC.csv', mode='a', sep=',', header=False, index=False)

#drop all models where >1 model is preferred
best = pd.read_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/best_AIC.csv')
dup_drop = best.drop_duplicates(subset= 'pair', keep=False)
dup_drop.to_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/best_AIC_dupdrop.csv', sep=',')

#keep only models where >1 model is preferred
dups = best[best.duplicated(subset='pair', keep = False)]
dups.to_csv('/Users/kericlamb/Desktop/Documents/Research/Spring_2021/DADI/All_results/res_vcf/best_AIC_duponly.csv', sep=',')



