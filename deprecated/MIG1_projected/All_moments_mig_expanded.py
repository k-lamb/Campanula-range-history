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

PMmod=open('./outputs/%s_MIG_expanded_output.txt' % Pair_name,'w')
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
    str("Tsc")+'\t'+ #Tsc - NA here
    str("m12")+'\t'+ #m12
    str("m21")+'\t'+ #m21
    str("theta")+'\t'+
    str("ll_model")+'\t'+
    str("aic")+'\n')
PMmod.close()

#model with symmetric migration from Moments bitbucket
def IM(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    nu1, nu2, Ts, m12, m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], Ts, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.pop_ids = pop_ids
    return fs

# split where one population can grow. from Momigliano et al. 2021
def IM_b(params, ns):
    nu1,nu2,s,Ts,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ts)
    nu_func= lambda t: [nu1,nu2_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, Ts, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    return fs

#reverse model to account for which pop expands
def IM_br(params, ns):
    nu1,nu2,s,Ts,m12,m21 = params
    nu1_0 = nu2*s
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ts)
    nu_func= lambda t: [nu2,nu1_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, Ts, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    return fs


# ancestral population has instant growth event 
def IM_ae(params, ns):
    nu1,nu2,nu_ae,Ts,Tae,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], Ts, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    return fs

#instant growth of ancestral and growth of descendant population 
def IM_ae_b(params, ns):
    nu1,nu2,nu_ae,s,Ts,Tae,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ts)
    nu_func= lambda t: [nu1,nu2_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, Ts, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    return fs

#reverse model to account for which population grows
def IM_ae_br(params, ns):
    nu1,nu2,nu_ae,s,Ts,Tae,m12,m21 = params
    nu1_0 = nu2*s
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ts)
    nu_func= lambda t: [nu2,nu1_func(t)]
    # calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, Ts, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    return fs

#running the rats nest of models

#parameters

# IM
param_IM = ["nu1", "nu2", "Ts", "m12", "m21"]
upper_bound_IM = [20, 20, 10, 20, 20]
lower_bound_IM = [1e-3, 1e-3, 1e-3, 1e-3, 1e-3]

# IM_b
param_IMb = ["nu1", "nu2", "s", "Ts", "m12", "m21"]
upper_bound_IMb = [20,20,0.999,10,20,20]
lower_bound_IMb = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

#IM_br
param_IMbr = ["nu1", "nu2", "s", "Ts", "m12", "m21"]
upper_bound_IMbr = [20,20,0.999,10,20,20]
lower_bound_IMbr = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

# IM_ae
param_IMae = ["nu1", "nu2", "nu_ae", "Ts", "Tae", "m12", "m21"]
upper_bound_IMae = [20,20,20,10,10,20,20]
lower_bound_IMae = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

# IM_ae_b
param_IMaeb = ["nu1", "nu2", "nu_ae", "s", "Ts", "Tae", "m12", "m21"]
upper_bound_IMaeb = [20,20,20,0.999,10,10,20,20]
lower_bound_IMaeb = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

# IM_ae_br
param_IMaebr = ["nu1", "nu2", "nu_ae", "s", "Ts", "Tae", "m12", "m21"]
upper_bound_IMaebr = [20,20,20,0.999,10,10,20,20]
lower_bound_IMaebr = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]

param_list = [param_IM, param_IMb, param_IMbr, param_IMae, param_IMaeb, param_IMaebr]
upper_bound_list = [upper_bound_IM, upper_bound_IMb, upper_bound_IMbr, upper_bound_IMae, upper_bound_IMaeb, upper_bound_IMaebr]
lower_bound_list = [lower_bound_IM, lower_bound_IMb, lower_bound_IMbr, lower_bound_IMae, lower_bound_IMaeb, lower_bound_IMaebr]

for n in range(int(param_list)):
    func_moments = param_list[n]
    upper_bound = upper_bound_list[n]
    lower_bound = lower_bound_list[n]
    model_name = str(func_moments)

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
        
        if model_name == "IM":
            nu_ae = "NA"
            s = "NA"
            Ts = popt[2]
            Tae = "NA"
            m12 = popt[3]
            m21 = popt[4]
            rev="NA"
        elif model_name == "IM_b":
            nu_ae = "NA"
            s = popt[2]
            Ts = popt[3]
            Tae = "NA"
            m12 = popt[4]
            m21 = popt[5]
            rev="FALSE"
        elif model_name == "IM_br":
            nu_ae = "NA"
            s = popt[2]
            Ts = popt[3]
            Tae = "NA"
            m12 = popt[4]
            m21 = popt[5]
            rev = "TRUE"
        elif model_name == "IM_ae":
            nu_ae = popt[2]
            s = "NA"
            Ts = popt[3]
            Tae =popt[4]
            m12 = popt[5]
            m21 = popt[6]
            rev="NA"
        elif model_name == "IM_ae_b":
            nu_ae = popt[2]
            s = popt[3]
            Ts = popt[4]
            Tae =popt[5]
            m12 = popt[6]
            m21 = popt[7]
            rev="FALSE"
        elif model_name == "IM_ae_br":
            nu_ae = popt[2]
            s = popt[3]
            Ts = popt[4]
            Tae =popt[5]
            m12 = popt[6]
            m21 = popt[7]
            rev="TRUE"

        PMmod=open('./outputs/%s_MIG_expanded_output.txt' % Pair_name,'a')
        PMmod.write(
            str(Pair_name)+'\t'+ #pair name
            str(rev)+'\t'+ #whether it is reversed or not
            str(model_name)+'\t'+ #model name
            str(nu1)+'\t'+ #nu1
            str(nu2)+'\t'+ #nu2
            str(nu_ae)+'\t'+ #nu_ae
            str(s)+'\t'+ #s
            str(Ts)+'\t'+ #Ts
            str(Tae)+'\t'+ #Tae
            str("NA")+'\t'+ #Tsc - NA here
            str(m12)+'\t'+ #m12
            str(m21)+'\t'+ #m21
            str(theta)+'\t'+
            str(ll_model)+'\t'+
            str(aic)+'\n')
        PMmod.close()
    print("Moments finished running %s" % model_name)
print("Moments finished running")
