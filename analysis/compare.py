#! /home/alan/opt/anaconda3/bin/python3.7

import numpy as np

cluster_mag = np.loadtxt("../cluster_method/magz_l_32_beta_0.200000_t_40.txt")
markovc_mag = np.loadtxt("../mc_sampling/magz_l_32_beta_0.200000_t_40.txt")


cl_mean = np.mean(cluster_mag[:,0])
cl_err = np.std(cluster_mag[:,0])/np.sqrt(np.shape(cluster_mag)[0])

print("T=%d %.6e +- %.6e"%(0,cl_mean,cl_err))


nT = np.shape(markovc_mag)[1]
for t in range(nT):
    mk_mean = np.mean(markovc_mag[:,t])
    mk_err = np.std(markovc_mag[:,t])/np.sqrt(np.shape(markovc_mag)[0])

    cl_mean = np.mean(cluster_mag[:,t+1])
    cl_err = np.std(cluster_mag[:,t+1])/np.sqrt(np.shape(cluster_mag)[0])

    print("T=%d %.6e +- %.6e         %.6e +- %.6e "%(t,cl_mean,cl_err,mk_mean,mk_err))


