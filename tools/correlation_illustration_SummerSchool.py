#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Illustration of correlation in multivariate 
@Author: Janus Rønn Lind Kobbersmed, janus@cfin.au.dk or januslind@gmail.com
"""


import numpy as np
from matplotlib import pyplot as plt

figpath = '\\Users\\Lara Godbersen\\Documents\\GitHub\\Masters-thesis\\Plots\\Visualize connectivity\\'

n_T = 50 #Number of time points

ylim = 2.6 #Limit of y-axis in plot

sim_mean = np.array([0,0]) #Means of the two examples

cov = 1 #Choice of cov for illustation. Equal to correlation since var=1.

covmat_negcor = np.array([[1, -cov], [-cov, 1]])
covmat_poscor = np.array([[1, cov], [cov, 1]])
covmat_nocor = np.array([[1, 0], [0, 1]])


##Generate random sample from multivariate normal distribution
rng = np.random.default_rng(seed=123)

ts_nocor = rng.multivariate_normal(sim_mean, covmat_nocor, size = n_T)
ts_negcor = rng.multivariate_normal(sim_mean, covmat_negcor, size = n_T)
ts_poscor = rng.multivariate_normal(sim_mean, covmat_poscor, size = n_T)


###Plots
figname = 'ts_nocor.png'
fig, ax = plt.subplots(2,1)
for i in [0,1]:
    ax[i].plot(np.arange(n_T), ts_nocor[:,i])
    ax[i].set_ylim([-ylim,ylim])
plt.savefig(figpath + figname)
plt.close('all')

figname = 'ts_negcor.png'
fig, ax = plt.subplots(2,1)
for i in [0,1]:
    ax[i].plot(np.arange(n_T), ts_negcor[:,i])
    ax[i].set_ylim([-ylim,ylim])
plt.savefig(figpath + figname)
plt.close('all')

figname = 'ts_poscor.png'
fig, ax = plt.subplots(2,1)
for i in [0,1]:
    ax[i].plot(np.arange(n_T), ts_poscor[:,i])
    ax[i].set_ylim([-ylim,ylim])
plt.savefig(figpath + figname)
plt.close('all')