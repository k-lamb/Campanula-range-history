import os
from os import listdir
from os.path import isfile, join
import numpy as np
import pylab
import sys
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
#setting current working directory

'''
needed modifications:
 - concatenate all no_mig models into one sheet, etc. 
 - unsure whether to include headers or not. only include if possible to control spacing
'''

#Model aic generator TRANSLATED
model_caps = sys.argv[1] #i.e. "NO_MIG", "MIG1", "SC1", "AM", "ONE_POP"

#constants
mu=2.8e-9
L=8.3e3

#grabs all files in directory of interest and makes it a .txt df
#for one-pop outputs only
outs = [f for f in listdir('/home/ksl2za/automation_dadi/All/{0}_projected/outputs'.format(model_caps)) if isfile(join('/home/ksl2za/automation_dadi/All/{0}_projected/outputs'.format(model_caps), f))]
output_list = pd.DataFrame(outs)
output_list.to_csv('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_outputs_list.txt'.format(model_caps), index=False, header=False)

#imports list of output files to use from above
out_list = []
with open ('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_outputs_list.txt'.format(model_caps), 'rt') as out_file:
    for myline in out_file:               
        out_list.append(myline.rstrip(""'\n').strip('""').strip('txt').strip('.'))

if model_caps == "NO_MIG":
    PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'w')
    PMmod.write(
        str("index")+' '+
        str("pair")+' '+
        str("nu1")+' '+
        str("nu2")+' '+
        str("Ts")+' '+
        str("theta")+' '+
        str("ll_model")+' '+
        str("aic")+' '+
        str("nu1_num")+' '+
        str("nu2_num")+' '+
        str("dv_time")
        )
    PMmod.close()
            
    #imports actual output files. uses for loop to grab lowest aic run from each output file
    for i in range(len(out_list)):
        df = pd.read_csv('/home/ksl2za/automation_dadi/All/{0}_projected/outputs/{1}.txt'.format(model_caps, out_list[i]), sep="\t")
        df['pop1'] =  (df['nu1'] * (df['theta']/(4*mu*L)))
        df['pop2'] =  (df['nu2'] * (df['theta']/(4*mu*L)))
        df['div_time'] = (2*(df['theta']/(4*mu*L))*df['Ts'])
        aic = df[df.aic==(df.aic.min())] #used in lieu of <2 as all models have nearly same score
        #aic = df[df.aic<(df.aic.min()+2)] #grabs rows with lowest aic (<2 apart)
        df_dict = dict.fromkeys(aic.columns, '') 
        aic = aic.rename(columns = df_dict) #sets column headers to '' so it doesn't import them
        if df.empty:
            continue
        PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'a') #writes lowest aic score to file and appends for each run
        PMmod.write(
            str(aic)+'\n')
        PMmod.close()

elif model_caps == "MIG1":
    PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'w')
    PMmod.write(
        str("index")+' '+
        str("pair")+' '+
        str("nu1")+' '+
        str("nu2")+' '+
        str("Ts")+' '+
        str("m12")+' '+
        str("m21")+' '+
        str("theta")+' '+
        str("ll_model")+' '+
        str("aic")+' '+
        str("nu1_num")+' '+
        str("nu2_num")+' '+
        str("mji")+' '+
        str("mij")+' '+
        str("mig_pop1")+' '+
        str("mig_pop2")+' '+
        str("dv_time")
        )
    PMmod.close()
            
    #imports actual output files. uses for loop to grab lowest aic run from each output file
    for i in range(len(out_list)):
        df = pd.read_csv('/home/ksl2za/automation_dadi/All/{0}_projected/outputs/{1}.txt'.format(model_caps, out_list[i]), sep="\t")
        df['pop1'] =  (df['nu1'] * (df['theta']/(4*mu*L)))
        df['pop2'] =  (df['nu2'] * (df['theta']/(4*mu*L)))
        df['mij'] = (df['m12']/(2*(df['theta']/(4*mu*L))))
        df['mji'] = (df['m21']/(2*(df['theta']/(4*mu*L))))
        df['mig_pop1'] = df['mij']*(df['nu1']*(df['theta']/(4*mu*L)))
        df['mig_pop2'] = df['mji']*(df['nu2']*(df['theta']/(4*mu*L)))
        df['div_time'] = (2*(df['theta']/(4*mu*L))*df['Ts'])
        aic = df[df.aic==(df.aic.min())] #used in lieu of <2 as all models have nearly same score
        #aic = df[df.aic<(df.aic.min()+2)] #grabs rows with lowest aic (<2 apart)
        df_dict = dict.fromkeys(aic.columns, '') 
        aic = aic.rename(columns = df_dict) #sets column headers to '' so it doesn't import them
        if df.empty:
            continue
        PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'a') #writes lowest aic score to file and appends for each run
        PMmod.write(
            str(aic)+'\n')
        PMmod.close()

elif model_caps == "SC1":
    PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'w')
    PMmod.write(
        str("index")+' '+
        str("pair")+' '+
        str("nu1")+' '+
        str("nu2")+' '+
        str("Ts")+' '+
        str("Tsc")+' '+
        str("m12")+' '+
        str("m21")+' '+
        str("theta")+' '+
        str("ll_model")+' '+
        str("aic")+' '+
        str("nu1_num")+' '+
        str("nu2_num")+' '+
        str("mji")+' '+
        str("mij")+' '+
        str("mig_pop1")+' '+
        str("mig_pop2")+' '+
        str("dv_time")+' '+
        str("sc_time")
        )
    PMmod.close()
            
    #imports actual output files. uses for loop to grab lowest aic run from each output file
    for i in range(len(out_list)):
        df = pd.read_csv('/home/ksl2za/automation_dadi/All/{0}_projected/outputs/{1}.txt'.format(model_caps, out_list[i]), sep="\t")
        df['pop1'] =  (df['nu1'] * (df['theta']/(4*mu*L)))
        df['pop2'] =  (df['nu2'] * (df['theta']/(4*mu*L)))
        df['mij'] = (df['m12']/(2*(df['theta']/(4*mu*L))))
        df['mji'] = (df['m21']/(2*(df['theta']/(4*mu*L))))
        df['mig_pop1'] = df['mij']*(df['nu1']*(df['theta']/(4*mu*L)))
        df['mig_pop2'] = df['mji']*(df['nu2']*(df['theta']/(4*mu*L)))
        df['div_time'] = (2*(df['theta']/(4*mu*L))*df['Ts'])
        df['sc_time'] = (2*(df['theta']/(4*mu*L))*df['Tsc'])
        aic = df[df.aic==(df.aic.min())] #used in lieu of <2 as all models have nearly same score
        #aic = df[df.aic<(df.aic.min()+2)] #grabs rows with lowest aic (<2 apart)
        df_dict = dict.fromkeys(aic.columns, '') 
        aic = aic.rename(columns = df_dict) #sets column headers to '' so it doesn't import them
        if df.empty:
            continue
        PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'a') #writes lowest aic score to file and appends for each run
        PMmod.write(
            str(aic)+'\n')
        PMmod.close()

elif model_caps == "AM":
    PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'w')
    PMmod.write(
        str("index")+' '+
        str("pair")+' '+
        str("nu1")+' '+
        str("nu2")+' '+
        str("Ts")+' '+
        str("Tm")+' '+
        str("m12")+' '+
        str("m21")+' '+
        str("theta")+' '+
        str("ll_model")+' '+
        str("aic")+' '+
        str("nu1_num")+' '+
        str("nu2_num")+' '+
        str("mji")+' '+
        str("mij")+' '+
        str("mig_pop1")+' '+
        str("mig_pop2")+' '+
        str("dv_time")+' '+
        str("sc_time")
        )
    PMmod.close()
            
    #imports actual output files. uses for loop to grab lowest aic run from each output file
    for i in range(len(out_list)):
        df = pd.read_csv('/home/ksl2za/automation_dadi/All/{0}_projected/outputs/{1}.txt'.format(model_caps, out_list[i]), sep="\t")
        df['pop1'] =  (df['nu1'] * (df['theta']/(4*mu*L)))
        df['pop2'] =  (df['nu2'] * (df['theta']/(4*mu*L)))
        df['mij'] = (df['m12']/(2*(df['theta']/(4*mu*L))))
        df['mji'] = (df['m21']/(2*(df['theta']/(4*mu*L))))
        df['mig_pop1'] = df['mij']*(df['nu1']*(df['theta']/(4*mu*L)))
        df['mig_pop2'] = df['mji']*(df['nu2']*(df['theta']/(4*mu*L)))
        df['div_time'] = (2*(df['theta']/(4*mu*L))*df['Ts'])
        df['am_time'] = (2*(df['theta']/(4*mu*L))*df['Tm'])
        aic = df[df.aic==(df.aic.min())] #used in lieu of <2 as all models have nearly same score
        #aic = df[df.aic<(df.aic.min()+2)] #grabs rows with lowest aic (<2 apart)
        df_dict = dict.fromkeys(aic.columns, '') 
        aic = aic.rename(columns = df_dict) #sets column headers to '' so it doesn't import them
        if df.empty:
            continue
        PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'a') #writes lowest aic score to file and appends for each run
        PMmod.write(
            str(aic)+'\n')
        PMmod.close()

else: #must be model_caps == ONE_POP            
    PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'w')
    PMmod.write(
        str("index")+' '+
        str("pair")+' '+
        str("theta")+' '+
        str("ll_model")+' '+
        str("aic")+'\n')
    PMmod.close()
            
    #imports actual output files. uses for loop to grab lowest aic run from each output file
    for i in range(len(out_list)):
        df = pd.read_csv('/home/ksl2za/automation_dadi/All/ONE_POP_projected/outputs/%s.txt' % out_list[i], sep="\t")
        aic = df[df.index==(df.index.min())] #used in lieu of aic as all models have same aic score
        df_dict = dict.fromkeys(aic.columns, '') 
        aic = aic.rename(columns = df_dict) #sets column headers to '' so it doesn't import them
        if df.empty:
            continue
        PMmod=open('/home/ksl2za/automation_dadi/All/ONE_POP_projected/ONE_POP_aic_max_best_translated.txt','a') #writes lowest aic score to file and appends for each run
        PMmod.write(
            str(aic)+'\n')
        PMmod.close()

    #imports actual output files. uses for loop to grab lowest aic run from each output file
    for i in range(len(out_list)):
        headerList = ['pair', 'theta', 'll_model', 'aic']
        df = pd.read_csv('/home/ksl2za/automation_dadi/All/{0}_projected/outputs/{1}.txt'.format(model_caps, out_list[i]), names=headerList, sep="\t")
        aic = df[df.aic==(df.aic.min())] #used in lieu of <2 as all models have nearly same score
        #aic = df[df.aic<(df.aic.min()+2)] #grabs rows with lowest aic (<2 apart)
        df_dict = dict.fromkeys(aic.columns, '') 
        aic = aic.rename(columns = df_dict) #sets column headers to '' so it doesn't import them
        if df.empty:
            continue
        PMmod=open('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps),'a') #writes lowest aic score to file and appends for each run
        PMmod.write(
            str(aic)+'\n')
        PMmod.close()

#GRABBING best (aic) model of each pair for each demographic scenario and putting it into a file for TRANSLATED
pair_name = []
with open ('/home/ksl2za/automation_dadi/All/All_pairs.txt', 'rt') as pair_file:
    for myline in pair_file:               
        pair_name.append(myline.rstrip(""'\n'))

#put all for a given model in one file and sort
df=pd.read_csv('/home/ksl2za/automation_dadi/All/{0}_projected/{0}_aic_max_best_translated.txt'.format(model_caps), delim_whitespace=True, index_col=0)
df=df.sort_values('pair')
df=df.drop_duplicates() #deals with duplicates caused by identity in one pop models
df.to_csv('/home/ksl2za/automation_dadi/All/{0}_best_sorted.csv'.format(model_caps), sep=',')
