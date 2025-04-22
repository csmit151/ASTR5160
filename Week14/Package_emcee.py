#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
emcee


@author: csmit151
"""

# CS Import packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join
from tqdm import tqdm
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from tqdm import tqdm
import random
import os
import emcee


# CS Impot line.data
line_data = np.loadtxt("/d/scratch/ASTR5160/week13/line.data", unpack=True)

# CS Global variables
mean = []
var = []

# CS Calculate mean and variance
for i in range(0,10,1):
    mean.append(np.mean(line_data[i]))
    var.append(np.var(line_data[i], ddof=1))
    
# CS Adapt turotrial to fit a linear modela to line_data

x = np.arange(1,11,1) # CS our x values
yerr = np.sqrt(var) # CS Recall sig**2 = var
y = mean

A = np.vander(x,2)
C = np.diag(yerr**2)
ATA = np.dot(A.T, A / (yerr**2)[:, None])
cov = np.linalg.inv(ATA)
w = np.linalg.solve(ATA, np.dot(A.T, y / yerr**2))
print("Least-squares estimates:")
print("m = {0:.3f} ± {1:.3f}".format(w[0], np.sqrt(cov[0, 0])))
print("b = {0:.3f} ± {1:.3f}".format(w[1], np.sqrt(cov[1, 1])))


# CS Close to LikehoodFunc values so this is working

def log_likelihood(theta, x, y, yerr):
    m, b, log_f = theta
    model = m * x + b
    sigma2 = yerr**2 + model**2 * np.exp(2 * log_f)
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

# CS Choose some inital values that make sense for our line
m_true = 3
b_true = 5
f_true = 0.5

from scipy.optimize import minimize

np.random.seed(42)
nll = lambda *args: -log_likelihood(*args)
initial = np.array([m_true, b_true, np.log(f_true)]) + 0.1 * np.random.randn(3)
soln = minimize(nll, initial, args=(x, y, yerr))
m_ml, b_ml, log_f_ml = soln.x

print("Maximum likelihood estimates:")
print("m = {0:.3f}".format(m_ml))
print("b = {0:.3f}".format(b_ml))
print("f = {0:.3f}".format(np.exp(log_f_ml)))

# CS Make our prior bounds realistic for our line
def log_prior(theta):
    m, b, log_f = theta
    if 0 < m < 7.0 and 0.0 < b < 10.0 and -10.0 < log_f < 1.0:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)


import emcee

pos = soln.x + 1e-4 * np.random.randn(32, 3)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(x, y, yerr)
)
sampler.run_mcmc(pos, 5000, progress=True);


fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["m", "b", "log(f)"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");

tau = sampler.get_autocorr_time()
print("Burn in:",tau)

flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
print("flat list of samples shape:",flat_samples.shape)

import corner

fig = corner.corner(flat_samples, labels=labels)

from IPython.display import display, Math

for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    #txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    #txt = txt.format(mcmc[1], q[0], q[1], labels[i])
    print(labels[i], "=",mcmc[1], "+/-:",q[0],",", q[1])

plt.show()

"""
Results
Least-squares estimates:
m = 3.026 ± 0.173
b = 3.172 ± 1.112
Maximum likelihood estimates:
m = 3.026
b = 3.172
f = 0.000
Burn in: [34.58449576 35.02090667 42.95929745]
flat list of samples shape: (10432, 3)
m = 3.028824908437665 +/-: 0.17672583027047217 , 0.17472596997128464
b = 3.132634664006285 +/-: 1.1012804819455115 , 1.1477778354557682
log(f) = -6.7697336454477135 +/-: 2.157652488112161 , 2.3445043484756543


"""