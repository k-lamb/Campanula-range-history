import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
import dadi
from dadi import Misc
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#define sys args
iterations = int(sys.argv[1])
Pair_name = sys.argv[2]
pop_id1 = sys.argv[3]
pop_id2 = sys.argv[4]
proj_n1 = sys.argv[5]
proj_n2 = sys.argv[6]

#converting floats to integers
proj_n1= int(proj_n1)
proj_n2= int(proj_n2)

#setting pop id's and projections
pop_id=[pop_id1,pop_id2]
ns=[proj_n1, proj_n2]

#read in data as file 
dd = dadi.Misc.make_data_dict_vcf("/scratch/ksl2za/VCFs/vcf_dp_noTperf_noVA112_missing725.recode.vcf", "/home/ksl2za/automation_dadi/All/batch_90_p25_NoCpMt_popmap.txt")
print("done reading in VCF")

fs = dadi.Spectrum.from_data_dict(dd, pop_ids= pop_id, projections = ns, polarized = False)
print("done formatting FS")

PMmod=open('./outputs/%s_ONEPOP_output.txt' % Pair_name,'w')
PMmod.write(
    str("Pair_name")+'\t'+ #pair name
    str("theta")+'\t'+
    str("ll_model")+'\t'+
    str("aic")+'\n')
PMmod.close()

#snm model from moments bitbucket
def snm(params, pop_ids=None):

    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.pop_ids = pop_ids
    return fs


func_moments = snm

for i in range(int(iterations)):
    print("starting optimization "+str(i))
    popt=[]
    params = 0
    model=func_moments(popt, ns)
    ll_model=moments.Inference.ll_multinom(model, fs)
    aic = 2*params - 2*ll_model
    print('Maximum log composite likelihood: {0}'.format(ll_model))
    theta = moments.Inference.optimal_sfs_scaling(model, fs)

    PMmod=open('./outputs/%s_ONEPOP_output.txt' % Pair_name,'a')
    PMmod.write(
    	str(Pair_name)+'\t'+ #pair name
        str(theta)+'\t'+
        str(ll_model)+'\t'+
        str(aic)+'\n')
    PMmod.close()
print("Moments finished running")
