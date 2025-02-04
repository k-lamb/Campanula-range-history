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
dd = dadi.Misc.make_data_dict_vcf("/scratch/ksl2za/demo_infer/65_linreg/noTperf_65miss_50miss_FULL.recode.vcf", "/scratch/ksl2za/demo_infer/65_linreg/vcf_popmap_linreg.txt")
print("done reading in VCF")

fs = dadi.Spectrum.from_data_dict(dd, pop_ids= pop_id, projections = ns, polarized = False)
print("done formatting FS")

PMmod=open('./outputs/%s_AM_expanded_output.txt' % Pair_name,'w')
PMmod.write(
    str("Pair_name")+'\t'+ #pair name
    str("reversed")+'\t'+ #whether pair is reversed in model or not
    str("model")+'\t'+ #model name
    str("nu1")+'\t'+ #nu1
    str("nu2")+'\t'+ #nu2
    str("nu_ae")+'\t'+ #nu_ae
    str("s")+'\t'+ #s
    str("Ts")+'\t'+ #Ts
    str("Tae")+'\t'+ #Tae
    str("Tsc")+'\t'+ #Tsc
    str("m12")+'\t'+ #m12
    str("m21")+'\t'+ #m21
    str("theta")+'\t'+
    str("ll_model")+'\t'+
    str("aic")+'\n')
PMmod.close()

#model with symmetric migration from Moments bitbucket
def AM(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    nu1, nu2, Ts, Tsc, m12, m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], Ts, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1, nu2], Tsc, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.pop_ids = pop_ids
    return fs

def AM_b(params, ns):
    nu1,nu2,s,Ts,Tsc,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Tsc)
    nu_func= lambda t: [nu1,nu2_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], Ts, dt_fac=0.01, m = np.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, Tsc, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    return fs

def AM_br(params, ns):
    nu1,nu2,s,Ts,Tsc,m12,m21 = params
    nu1_0 = nu2*s
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Tsc)
    nu_func= lambda t: [nu2,nu1_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], Ts, dt_fac=0.01, m = np.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, Tsc, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    return fs

def AM_ae(params, ns):
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], Ts, dt_fac=0.01, m = np.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1, nu2], Tsc, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    return fs

def AM_ae_b(params, ns):
    nu1,nu2,nu_ae,s,Tae,Ts,Tsc,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Tsc)
    nu_func= lambda t: [nu1,nu2_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], Ts, dt_fac=0.01, m = np.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, Tsc, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    return fs

def AM_ae_br(params, ns):
    nu1,nu2,nu_ae,s,Tae,Ts,Tsc,m12,m21 = params
    nu1_0 = nu2*s
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Tsc)
    nu_func= lambda t: [nu2,nu1_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], Ts, dt_fac=0.01, m = np.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, Tsc, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    return fs

#parameters 

#AM
param_AM = ["nu1", "nu2", "Ts", "Tsc", "m12", "m21"]
upper_bound_AM = [10, 10, 5, 5, 20, 20]
lower_bound_AM = [1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3]

#AM_b
param_AMb = ["nu1", "nu2", "s", "Ts", "Tsc", "m12", "m21"]
upper_bound_AMb = [20,20,0.999,10,10,20,20]
lower_bound_AMb = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

#AM_ae
param_AMae = ["nu1", "nu2", "nu_ae", "Tae", "Ts", "Tsc", "m12", "m21"]
upper_bound_AMae = [20,20,20,10,10,10,20,20]
lower_bound_AMae = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

#AM_ae_b
param_AMaeb = ["nu1", "nu2", "nu_ae", "s", "Tae", "Ts", "Tsc", "m12", "m21"]
upper_bound_AMaeb = [20,20,20,0.999,10,10,10,20,20]
lower_bound_AMaeb = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

param_list = [param_AM, param_AMb, param_AMb, param_AMae, param_AMaeb, param_AMaeb]
upper_bound_list = [upper_bound_AM, upper_bound_AMb, upper_bound_AMb, upper_bound_AMae, upper_bound_AMaeb, upper_bound_AMaeb]
lower_bound_list = [lower_bound_AM, lower_bound_AMb, lower_bound_AMb, lower_bound_AMae, lower_bound_AMaeb, lower_bound_AMaeb]
model_list = [AM, AM_b, AM_br, AM_ae, AM_ae_b, AM_ae_br]

for n in range(len(param_list)):
    func_moments = model_list[n]
    upper_bound = upper_bound_list[n]
    lower_bound = lower_bound_list[n]
    model_name = model_list[n]

    for i in range(int(iterations)):
        print("starting optimization "+str(i)+"for %s" % model_name)
        params = len(param_list[n])
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
            model_name_abbrev = "AM"
            nu_ae = "NA"
            s = "NA"
            Tae = "NA"
            Ts = popt[2]
            Tsc = popt[3]
            m12 = popt[4]
            m21 = popt[5]
            rev = "NA"
        elif model_name == model_list[1]:
            model_name_abbrev = "AM_b"
            nu_ae = "NA"
            s = popt[2]
            Tae = "NA"
            Ts = popt[3]
            Tsc = popt[4]
            m12 = popt[5]
            m21 = popt[6]
            rev = "FALSE"
        elif model_name == model_list[2]:
            model_name_abbrev = "AM_br"
            nu_ae = "NA"
            s = popt[2]
            Tae = "NA"
            Ts = popt[3]
            Tsc = popt[4]
            m12 = popt[5]
            m21 = popt[6]
            rev = "TRUE"
        elif model_name == model_list[3]:
            model_name_abbrev = "AM_ae"
            nu_ae = popt[2]
            s = "NA"
            Tae = popt[3]
            Ts = popt[4]
            Tsc = popt[5]
            m12 = popt[6]
            m21 = popt[7]
            rev = "NA"
        elif model_name == model_list[4]:
            model_name_abbrev = "AM_ae_b"
            nu_ae = popt[2]
            s = popt[3]
            Tae = popt[4]
            Ts = popt[5]
            Tsc = popt[6]
            m12 = popt[7]
            m21 = popt[8]
            rev = "FALSE"
        elif model_name == model_list[5]:
            model_name_abbrev = "AM_ae_br"
            nu_ae = popt[2]
            s = popt[3]
            Tae = popt[4]
            Ts = popt[5]
            Tsc = popt[6]
            m12 = popt[7]
            m21 = popt[8]
            rev = "TRUE"

        PMmod=open('./outputs/%s_AM_expanded_output.txt' % Pair_name,'a')
        PMmod.write(
            str(Pair_name)+'\t'+ #pair name
            str(rev)+'\t'+ #reversed pairnames, true or false
            str(model_name_abbrev)+'\t'+ #model name
            str(nu1)+'\t'+ #nu1
            str(nu2)+'\t'+ #nu2
            str(nu_ae)+'\t'+ #nu_ae
            str(s)+'\t'+ #s
            str(Ts)+'\t'+ #Ts <- time of migration and split 
            str(Tae)+'\t'+ #Tae 
            str(Tsc)+'\t'+ #Tsc <- time migration ceases 
            str(m12)+'\t'+ #m12
            str(m21)+'\t'+ #m21
            str(theta)+'\t'+
            str(ll_model)+'\t'+
            str(aic)+'\n')
        PMmod.close()
    print("Moments finished running %s" % model_name)
print("Moments finished running")
