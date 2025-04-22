#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Likelihood Functions and MCMC

@author: csmit151
"""

# CS Import packages
import numpy as np
import random


# CS Impot line.data
line_data = np.loadtxt("/d/scratch/ASTR5160/week13/line.data", unpack=True)

# CS Global variables
mean = []
var = []

# CS Calculate mean and variance
for i in range(0,10,1):
    mean.append(np.mean(line_data[i]))
    var.append(np.var(line_data[i], ddof=1))

# CS Write a funciton to calculate (ln) posterior probability
# CS for a linear fit to the data when passed m and b

def LogLikeLine(m,b,data):
    """ Takes m, b, and data for y=mx+b
    and calculates posterior probability    
    
    Parameters
    ----------
    m: :class: 'float'
        Slope value
    b: :class: 'float'
        Intercept value
    data: :class: 'array'
        Data of drawn from a line of the form y=mx+b    
        
    Returns
    -------
    Posterior probability value for a given m and b    
    """
    to_sum = []
    mean = []
    var = []
    for i in range(len(data)): # CS Calculate mean and variance
        mean.append(np.mean(data[i]))
        var.append(np.var(data[i], ddof=1))
    for k in range(len(data)): # CS Calculate posterior probability, assuming x value of data corresponds to 0<x<1 (bin x_0)
        to_sum.append(((mean[k]-(m*(k+1)+b))**2)/var[k]+np.log(2*np.pi*var[k]**2))
    
    # CS Calculate logL
    logL = -0.5 * np.sum(to_sum)
    # CS Apply our flat prior
    if (0 < b < 8):
        prior = 0
    else:
        prior = -np.inf # CS ln(0) = -inf
    return logL + prior


def MHwalk(step_size,m_0,b_0,data,steps):
    """ Does a MH walk through a linear parameter space
    and returns an MCMC chain and its best fit values
    
    Parameters
    ----------
    step_zise: :class: 'float'
        The standard deviation given to the Gaussian proposal function
    m_0: :class: 'float'
        Starting guess of m
    b_0: :class: 'float'
        Starting guess of b
     data: :class: 'array'
         Data of drawn from a line of the form y=mx+b
     steps: :class: 'float'
         Number of MH steps (i.e. length of MCMC chain)   
    
    Returns
    -------
    Best fit m, b, and maximum chain value
    """
    # CS Calculate starting parameters
    chain = []
    mchain = []
    bchain = []
    for i in range(0,steps,1):
        LL = LogLikeLine(m_0, b_0, data)
        del_m = m_0-random.gauss(m_0,step_size) # CS A Gaussian centered on m_0, where we get a random value of del_m
        del_b = m_0-random.gauss(b_0,step_size)
        m_new = m_0+del_m
        b_new = b_0+del_b
        LL1 = LogLikeLine(m_new, b_new, data)
        # CS Calculate R
        R = np.e**(LL1-LL)
        if R > 1:
            # CS Accept the new parameters
            chain.append(LL1)
            mchain.append(m_new)
            bchain.append(b_new)
            m_0 = m_new
            b_0 = b_new
        if R < 1:
            if np.random.rand() < R: # Accept new parameters with probability R
                chain.append(LL1)
                mchain.append(m_new)
                bchain.append(b_new)
                m_0 = m_new
                b_0 = b_new
    ii = np.argmax(chain) # CS Chain values corresponding to the max pprob represent model that best fits the data
    print("Best fit m:",mchain[ii])
    print("Best fit b:",bchain[ii])
    print("Best fit Max pp:",chain[ii])


if __name__ == "__main__":
    ss = 0.1
    mstart = 3
    bstart = 4
    Nsteps = 200
    MHwalk(ss,mstart,bstart,line_data,Nsteps)




