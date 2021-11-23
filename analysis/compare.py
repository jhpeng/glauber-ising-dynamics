#! /home/alan/opt/anaconda3/bin/python3.7

import numpy as np
import sys

L=int(sys.argv[1])
t=int(sys.argv[2])
mz=float(sys.argv[3])
beta=float(sys.argv[4])

cluster_mag = np.loadtxt("../cluster_method/magz_l_%d_mz_%.6f_beta_%.6f_t_%d.txt"%(L,mz,beta,t))
markovc_mag = np.loadtxt("../mc_sampling/magz_l_%d_mz_%.6f_beta_%.6f_t_%d.txt"%(L,mz,beta,t))


cl_mean = np.mean(cluster_mag[:,0])
cl_err = np.std(cluster_mag[:,0])/np.sqrt(np.shape(cluster_mag)[0])

print("T=%d %.6e +- %.6e"%(0,cl_mean,cl_err))



nT = np.shape(markovc_mag)[1]
data = np.zeros([nT+1,5])


data[0,0] = 0;
data[0,1] = cl_mean;
data[0,2] = cl_err;
data[0,3] = 0;
data[0,4] = 0;

for t in range(nT):
    mk_mean = np.mean(markovc_mag[:,t])
    mk_err = np.std(markovc_mag[:,t])/np.sqrt(np.shape(markovc_mag)[0])

    cl_mean = np.mean(cluster_mag[:,t+1])
    cl_err = np.std(cluster_mag[:,t+1])/np.sqrt(np.shape(cluster_mag)[0])

    print("T=%d %.6e +- %.6e         %.6e +- %.6e "%(t+1,cl_mean,cl_err,mk_mean,mk_err))

    data[t+1,0] = t+1;
    data[t+1,1] = cl_mean;
    data[t+1,2] = cl_err;
    data[t+1,3] = mk_mean;
    data[t+1,4] = mk_err;

np.savetxt("magz_l_%d_mz_%.6f_beta_%.6f_t_%d.txt"%(L,mz,beta,t),data)

