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

#maxiter setting
maxi = 100

#read in data as file 
dd = dadi.Misc.make_data_dict_vcf("/scratch/ksl2za/Moments/vcf_dp_noTperf_noVA112_missing725.recode.vcf", "/scratch/ksl2za/Moments/batch_90_p25_NoCpMt_popmap.txt")
print("done reading in VCF")

fs = dadi.Spectrum.from_data_dict(dd, pop_ids= pop_id, projections = ns, polarized = False)
print("done formatting FS")

PMmod=open('./outputs/%s_NOMIG_expanded_output.txt' % Pair_name,'w')
PMmod.write(
    str("Pair_name")+'\t'+ #pair name
    str("reversed")+'\t'+ #whether pair is reversed or not
    str("model")+'\t'+ #model name
    str("nu1")+'\t'+ #nu1
    str("nu2")+'\t'+ #nu2
    str("nu_ae")+'\t'+ #nu_ae
    str("s")+'\t'+ #s
    str("Ts")+'\t'+ #Ts
    str("Tae")+'\t'+ #Tae
    str("Tsc")+'\t'+ #Tsc - NA here
    str("m12")+'\t'+ #m12 - NA here
    str("m21")+'\t'+ #m21 - NA here
    str("theta")+'\t'+
    str("ll_model")+'\t'+
    str("aic")+'\n')
PMmod.close()

#model with symmetric migration from Moments bitbucket
def NM(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    nu1, nu2, Ts = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], Ts, dt_fac=0.01)
    fs.pop_ids = pop_ids
    return fs

# split where one population can grow. from Momigliano et al. 2021
def NM_b(params, ns):
    nu1,nu2,s,Ts= params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ts)
    nu_func= lambda t: [nu1,nu2_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, Ts, dt_fac=0.01)
    return fs

# split where one population can grow. from Momigliano et al. 2021
def NM_br(params, ns):
    nu1,nu2,s,Ts= params
    nu1_0 = nu2*s
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ts)
    nu_func= lambda t: [nu2,nu1_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, Ts, dt_fac=0.01)
    return fs

# ancestral population has instant growth event 
def NM_ae(params, ns):
    nu1,nu2,nu_ae,Ts,Tae = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], Ts, dt_fac=0.01)
    return fs

#instant growth of ancestral and growth of descendant population 
def NM_ae_b(params, ns):
    nu1,nu2,nu_ae,s,Ts,Tae = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ts)
    nu_func= lambda t: [nu1,nu2_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, Ts, dt_fac=0.01)
    return fs

#instant growth of ancestral and growth of descendant population 
def NM_ae_b(params, ns):
    nu1,nu2,nu_ae,s,Ts,Tae = params
    nu1_0 = nu2*s
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ts)
    nu_func= lambda t: [nu2,nu1_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, Ts, dt_fac=0.01)
    return fs

#instant growth of ancestral and growth of descendant population 
def NM_ae_br(params, ns):
    nu1,nu2,nu_ae,s,Ts,Tae = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ts)
    nu_func= lambda t: [nu1,nu2_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, Ts, dt_fac=0.01)
    return fs

#parameters 

#NM
param_NM = ["nu1", "nu2", "Ts"]
upper_bound_NM = [10, 10, 5]
lower_bound_NM = [1e-3, 1e-3, 1e-3]

#NM_b
param_NMb = ["nu1", "nu2", "s", "Ts"]
upper_bound_NMb = [20,20,0.999,10]
lower_bound_NMb = [1e-3,1e-3,1e-3,1e-3]

#NM_ae
param_NMae = ["nu1", "nu2", "nu_ae", "Tae", "Ts"]
upper_bound_NMae = [20,20,20,10,10]
lower_bound_NMae = [1e-3,1e-3,1e-3,1e-3,1e-3]

#NM_ae_b
param_NMaeb = ["nu1", "nu2", "nu_ae", "s", "Tae", "Ts"]
upper_bound_NMaeb = [20,20,20,0.999,10,10]
lower_bound_NMaeb = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

param_list = [param_NM, param_NMb, param_NMb, param_NMae, param_NMaeb, param_NMaeb]
upper_bound_list = [upper_bound_NM, upper_bound_NMb, upper_bound_NMb, upper_bound_NMae, upper_bound_NMaeb, upper_bound_NMaeb]
lower_bound_list = [lower_bound_NM, lower_bound_NMb, lower_bound_NMb, lower_bound_NMae, lower_bound_NMaeb, lower_bound_NMaeb]
model_list = [NM, NM_b, NM_br, NM_ae, NM_ae_b, NM_ae_br]

for n in range(int(param_list)):
    func_moments = model_list[n]
    upper_bound = upper_bound_list[n]
    lower_bound = lower_bound_list[n]
    model_name = model_list[n]

    for i in range(len(iterations)):
        print("starting optimization "+str(i)+"for %s" % model_name)
        params = len(["nu1", "nu2", "Ts"])
        popt=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(params)]
        popt=moments.Inference.optimize_log(popt, fs, func_moments,
                                            lower_bound=lower_bound, upper_bound=upper_bound,
                                            verbose=False, maxiter=maxi,)
        model=func_moments(popt, ns)
        ll_model=moments.Inference.ll_multinom(model, fs)
        aic = 2*params - 2*ll_model
        print('Maximum log composite likelihood: {0}'.format(ll_model))
        theta = moments.Inference.optimal_sfs_scaling(model, fs)
        
        #sorting out issue with parameter recording in outputs
        nu1 = popt[0]
        nu2 = popt[1]
        
        if model_name == model_list[0]:
            model_name_abbrev = "NM"
            nu_ae = "NA"
            s = "NA"
            Tae = "NA"
            Ts = popt[2]
            rev = "NA"
        elif model_name == model_list[1]:
            model_name_abbrev = "NM_b"
            nu_ae = "NA"
            s = popt[2]
            Tae = "NA"
            Ts = popt[3]
            rev = "FALSE"
        elif model_name == model_list[2]:
            model_name_abbrev = "NM_br"
            nu_ae = "NA"
            s = popt[2]
            Tae = "NA"
            Ts = popt[3]
            rev = "TRUE"
        elif model_name == model_list[3]:
            model_name_abbrev = "NM_ae"
            nu_ae = popt[2]
            s = "NA"
            Tae = popt[3]
            Ts = popt[4]
            rev = "NA"
        elif model_name == model_list[4]:
            model_name_abbrev = "NM_ae_b"
            nu_ae = popt[2]
            s = popt[3]
            Tae = popt[4]
            Ts = popt[5]
            rev = "FALSE"
        elif model_name == model_list[5]:
            model_name_abbrev = "NM_ae_br"
            nu_ae = popt[2]
            s = popt[3]
            Tae = popt[4]
            Ts = popt[5]
            rev = "TRUE"

        PMmod=open('./outputs/%s_NOMIG_expanded_output.txt' % Pair_name,'a')
        PMmod.write(
            str(Pair_name)+'\t'+ #pair name
            str(rev)+'\t'+ #whether pair is reversed or not. NA/T/F
            str(model_name_abbrev)+'\t'+ #model name
            str(nu1)+'\t'+ #nu1
            str(nu2)+'\t'+ #nu2
            str(nu_ae)+'\t'+ #nu_ae
            str(s)+'\t'+ #s
            str(Ts)+'\t'+ #Ts
            str(Tae)+'\t'+ #Tae
            str("NA")+'\t'+ #Tsc 
            str("NA")+'\t'+ #m12
            str("NA")+'\t'+ #m21
            str(theta)+'\t'+
            str(ll_model)+'\t'+
            str(aic)+'\n')
        PMmod.close()
    print("Moments finished running %s" % model_name)
print("Moments finished running")
