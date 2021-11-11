#! /home/alan/opt/anaconda3/bin/python3.7

import numpy as np

cluster_mag = np.loadtxt("../cluster_method/magz_l_16_beta_0.250000_t_8.txt")
markovc_mag = np.loadtxt("../mc_sampling/magz_l_16_beta_0.250000_t_8.txt")



nblock = np.shape(cluster_mag)[0]
nT = np.shape(cluster_mag)[1]

mean = np.mean(cluster_mag,axis=0)
std = np.std(cluster_mag,axis=0)
for t in range(nT):
    print("T=%d %.6e +- %.6e"%(t,mean[t],std[t]/np.sqrt(nblock)))

print("-------------------------------------")


nblock = np.shape(markovc_mag)[0]
nT = np.shape(markovc_mag)[1]

mean = np.mean(markovc_mag,axis=0)
std = np.std(markovc_mag,axis=0)
for t in range(nT):
    print("T=%d %.6e +- %.6e"%(t+1,mean[t],std[t]/np.sqrt(nblock)))

