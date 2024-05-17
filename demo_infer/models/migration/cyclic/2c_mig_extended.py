from socketserver import ThreadingUnixDatagramServer
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
iterations = sys.argv[1]
iterations = int(iterations)+90
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
dd = dadi.Misc.make_data_dict_vcf("/scratch/ksl2za/demo_infer/80_linreg_2c/noTperf_80miss_50miss_FULL.recode.vcf", "/scratch/ksl2za/demo_infer/80_linreg_2c/vcf_popmap_linreg.txt")
print("done reading in VCF")

fs = dadi.Spectrum.from_data_dict(dd, pop_ids= pop_id, projections = ns, polarized = False)
print("done formatting FS")

# PMmod=open('./outputs/%s_MIG_expanded_output.txt' % Pair_name,'w')
# PMmod.write(
#     str("Pair_name")+'\t'+ #pair name
#     str("reversed")+'\t'+ #whether pair is reversed in model or not
#     str("model")+'\t'+ #model name
#     str("nu1")+'\t'+ #nu1
#     str("nu2")+'\t'+ #nu2
#     str("nu_ae")+'\t'+ #nu_ae
#     str("s")+'\t'+ #s
#     str("Ts")+'\t'+ #Ts
#     str("Tae")+'\t'+ #Tae
#     str("Tsc")+'\t'+ #Tsc - NA here
#     str("Ts2")+'\t'+ #Ts2
#     str("Tsc2")+'\t'+ #Tsc2 - NA here
#     str("m12")+'\t'+ #m12
#     str("m21")+'\t'+ #m21
#     str("m12_2")+'\t'+ #m12_2
#     str("m21_2")+'\t'+ #m21_2
#     str("theta")+'\t'+
#     str("ll_model")+'\t'+
#     str("aic")+'\n')
# PMmod.close()

#model with symmetric migration from Moments bitbucket
def IM_2C(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    nu1, nu2, Ts, Ts2, Tsc, m12, m21, m12_2, m21_2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], Ts, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1, nu2], Ts)
    fs.integrate([nu1, nu2], Tsc, dt_fac=0.01, m=np.array([[0, m12_2], [m21_2, 0]]))
    fs.pop_ids = pop_ids
    return fs

# ancestral population has instant growth event 
def IM_2C_ae(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    nu1, nu2, nu_ae, Tae, Ts, Ts2, Tsc, m12, m21, m12_2, m21_2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], Ts, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1, nu2], Ts)
    fs.integrate([nu1, nu2], Tsc, dt_fac=0.01, m=np.array([[0, m12_2], [m21_2, 0]]))
    fs.pop_ids = pop_ids
    return fs

#running the rats nest of models

#parameters

# IM
param_IM_2C = ["nu1", "nu2", "Ts", "Ts2", "Tsc", "m12", "m21", "m12_2", "m21_2"]
upper_bound_IM = [20, 20, 10, 10, 10, 20, 20, 20, 20]
lower_bound_IM = [1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3]

# IM_ae
param_IM_2C_ae = ["nu1", "nu2", "nu_ae", "Tae", "Ts", "Ts2", "Tsc", "m12", "m21", "m12_2", "m21_2"]
upper_bound_IMae = [20,20,20,10,10,10,10,20,20,20,20]
lower_bound_IMae = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

param_list = [param_IM_2C, param_IM_2C_ae]
upper_bound_list = [upper_bound_IM, upper_bound_IMae]
lower_bound_list = [lower_bound_IM, lower_bound_IMae]
model_list = [IM_2C, IM_2C_ae]

for n in range(len(param_list)):
    func_moments = model_list[n]
    upper_bound = upper_bound_list[n]
    lower_bound = lower_bound_list[n]
    model_name = model_list[n]

    for i in range(int(iterations)):
        #need to deal with the fact that we have to run IM_b and IM_ae_b twice to account for order of populations

        print("starting optimization "+str(i)+"for %s" % model_name)

        params = len(param_list[n])
        popt=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(params)]
        popt=moments.Inference.optimize_log(popt, fs, func_moments,
                                            lower_bound=lower_bound, upper_bound=upper_bound,
                                            verbose=False, maxiter=100,)
        model=func_moments(popt, ns)
        ll_model=moments.Inference.ll_multinom(model, fs)
        aic = 2*params - 2*ll_model
        print('Maximum log composite likelihood: {0}'.format(ll_model))
        theta = moments.Inference.optimal_sfs_scaling(model, fs)

        #sorting out issue with parameter recording in outputs
        nu1 = popt[0]
        nu2 = popt[1]
        if model_name == model_list[0]:
            model_name_abbrev = "IM_2C"
            nu_ae = "NA"
            s = "NA"
            Ts = popt[2]
            Tae = "NA"
            Ts2 = popt[3]
            Tsc = popt[4]
            Tsc2 = "NA"
            m12 = popt[5]
            m21 = popt[6]
            m12_2 = popt[7]
            m21_2 = popt[8]
            rev="NA"
        elif model_name == model_list[1]:
            model_name_abbrev = "IM_2C_ae"
            nu_ae = popt[2]
            s = "NA"
            Tae = popt[3]
            Ts = popt[4]
            Ts2 = popt[5]
            Tsc = popt[6]
            Tsc2 = "NA"
            m12 = popt[7]
            m21 = popt[8]
            m12_2 = popt[9]
            m21_2 = popt[10]
            rev="NA"

        PMmod=open('./outputs/%s_MIG_expanded_output.txt' % Pair_name,'a')
        PMmod.write(
            str(Pair_name)+'\t'+ #pair name
            str(rev)+'\t'+ #whether it is reversed or not
            str(model_name_abbrev)+'\t'+ #model name
            str(nu1)+'\t'+ #nu1
            str(nu2)+'\t'+ #nu2
            str(nu_ae)+'\t'+ #nu_ae
            str(s)+'\t'+ #s
            str(Ts)+'\t'+ #Ts
            str(Tae)+'\t'+ #Tae
            str(Tsc)+'\t'+ #Tsc
            str(Ts2)+'\t'+ #Ts2
            str(Tsc2)+'\t'+ #Tsc2
            str(m12)+'\t'+ #m12
            str(m21)+'\t'+ #m21
            str(m12_2)+'\t'+ #m12
            str(m21_2)+'\t'+ #m21
            str(theta)+'\t'+
            str(ll_model)+'\t'+
            str(aic)+'\n')
        PMmod.close()
    print("Moments finished running %s" % model_name)
print("Moments finished running")
